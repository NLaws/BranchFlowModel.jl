"""
    dsstxt_to_sparse_array(fp::String, first_data_row::Int = 5)

convert a SystemY.txt file from OpenDSS to a julia matrix.
assumes that Y is symmetric.
"""
function dsstxt_to_sparse_array(fp::String, first_data_row::Int = 5)

    rows = Int[]
    cols = Int[]
    real = Float64[]
    imag = Float64[]

    for (i, line) in enumerate(eachline(fp))

        if i < first_data_row continue end
        line = replace(line, " "=>"")  # "[1,1]=50500+j-50500"
        N = length(line)
        if N == 0 continue end

        append!(rows, tryparse(Int64,
                chop(line, head=findfirst("[", line)[end], tail=N-findfirst(",", line)[end]+1)
        ))

        append!(cols, tryparse(Int64,
                chop(line, head=findfirst(",", line)[end], tail=N-findfirst("]", line)[end]+1)
        ))

        append!(real, tryparse(Float64,
                chop(line, head=findfirst("=", line)[end], tail=N-findfirst("+", line)[end]+1)
        ))

        append!(imag, tryparse(Float64,
                chop(line, head=findfirst("j", line)[end], tail=0)
        ))
    end
    return convert(Array{Complex, 2}, Symmetric(sparse(rows, cols, complex.(real, imag)), :L))
end


function heads(edges:: Vector{Tuple})
    return collect(e[1] for e in edges)
end


function tails(edges:: Vector{Tuple})
    return collect(e[2] for e in edges)
end

"""
    dss_dict_to_arrays(d::Dict)

Parse the dict from PowerModelsDistribution.parse_dss into values needed for LinDistFlow
"""
function dss_dict_to_arrays(d::Dict)
    # TODO allocate empty arrays with number of lines
    # TODO separate this method into sub-methods, generally parse components separately, add transformers
    edges = Tuple[]
    phases = Vector[]
    linecodes = String[]
    linelengths = Float64[]
    Isqaured_up_bounds = Dict{String, Float64}()

    # some reuseable stuff
    function get_b1_b2_phs(v::Dict)
        if occursin(".", v["bus1"])  # have to account for .1.2 phases for example
            b1 = chop(v["bus1"], tail=length(v["bus1"])-findfirst('.', v["bus1"])+1)
            phs = sort!(collect(parse(Int,ph) for ph in split(v["bus1"][findfirst('.', v["bus1"])+1:end], ".")))
        else  # default to 3 phases
            b1 = v["bus1"]
            phs = [1,2,3]
        end

        if occursin(".", v["bus2"])
            b2 = chop(v["bus2"], tail=length(v["bus2"])-findfirst('.', v["bus2"])+1)
        else
            b2 = v["bus2"]
        end
        return b1, b2, phs
    end

    for (k,v) in d["line"]  # Line dict includes switches
        if "switch" in keys(v) && v["switch"] == true
            try
                # need to connect busses over switch
                b1, b2, phs = get_b1_b2_phs(v)
                linecode = "switch" * v["name"]
                push!(edges, (b1, b2))
                push!(linecodes, linecode)
                push!(linelengths, 1.0)
                push!(phases, phs)

                Isqaured_up_bounds[linecode] = DEFAULT_AMP_LIMIT^2
                if "emergamps" in keys(v)  # assuming lowercase keys
                    Isqaured_up_bounds[linecode] = v["emergamps"]^2
                elseif "normamps" in keys(v)
                    Isqaured_up_bounds[linecode] = v["normamps"]^2
                end

                # TODO handle missing r1 or x1
                d["linecode"][linecode] = Dict(
                    "nphases" => length(phs),
                    "rmatrix" => Diagonal(ones(3)) * v["r1"],
                    "xmatrix" => Diagonal(ones(3)) * v["x1"],
                )
            catch e
                @warn("Unable to parse switch $(k) when processing OpenDSS model.")
                println(e)
            end
            continue
        end
        try
            b1, b2, phs = get_b1_b2_phs(v)
            push!(edges, (b1, b2))
            push!(linecodes, v["linecode"])
            push!(phases, phs)

            # TODO ratings could be in linecode dict too
            # TODO assuming that there are linecodes, should converge on consistent keys for lines
            if "emergamps" in keys(v)  # assuming lowercase keys
                Isqaured_up_bounds[v["linecode"]] = v["emergamps"]^2
            elseif "normamps" in keys(v)
                Isqaured_up_bounds[v["linecode"]] = v["normamps"]^2
            else
                Isqaured_up_bounds[v["linecode"]] = DEFAULT_AMP_LIMIT^2
            end

            # TODO handle scaling of lengths and R/X values
            # for now just make sure the linecode and line values are in consistent units
            # and BEWARE PowerModelsDistribution will scale values from openDSS!
            push!(linelengths, v["length"]) 
        catch
            @warn("Unable to parse line $(k) when processing OpenDSS model.")
        end
    end

    # make phases_into_bus to infer transformer phases
    phases_into_bus = Dict(k=>v for (k,v) in zip(tails(edges), phases))

    for (k,v) in get(d, "transformer", Dict())
        try
            # need to connect busses over transformers
            b1 = get(v, "bus", nothing)
            b2 = get(v, "bus_2", nothing)

            if !(b1 in tails(edges)) && !(b2 in heads(edges))
                # this transformer does not connect anything so we ignore it
                continue
            end

            if b1 in tails(edges)
                phs = phases_into_bus[b1]
            else
                @warn("Not parsing transformer $k between $b1 and $b2
                      because it does not have a line or switch into it.")
                continue
            end

            nwindings = v["windings"]
            if nwindings != 2
                @warn("Parsing a $nwindings winding transformer as a 2 winding transformer.")
            end

            R1 = v["%r"] / 100 * v["kv"]^2 / v["kva"]
            R2 = v["%r_2"] / 100 * v["kv_2"]^2 / v["kva_2"]
            R = R1 + R2
            X = v["xhl"] / 100 * v["kv"]^2 / v["kva"]
            # TODO other reactance values XLT, XHT

            linecode = v["name"]
            push!(edges, (b1, b2))
            push!(linecodes, linecode)
            push!(linelengths, 1.0)
            push!(phases, phs)

            Isqaured_up_bounds[linecode] = DEFAULT_AMP_LIMIT^2

            rmatrix = zeros(3,3)
            xmatrix = zeros(3,3)
            # set the diagaonal values
            for phs1 in phs
                rmatrix[phs1, phs1] = R
                xmatrix[phs1, phs1] = X
            end

            # TODO handle missing r1 or x1
            d["linecode"][linecode] = Dict(
                "nphases" => length(phs),
                "rmatrix" => rmatrix,
                "xmatrix" => xmatrix,
            )
        catch e
            @warn("Unable to parse transformer $(k) when processing OpenDSS model.")
            println(e)
        end
    end

    return edges, linecodes, linelengths, d["linecode"], phases, Isqaured_up_bounds
end


"""
    dss_loads(d::Dict)

Return the P,Q loads from the dict provided by parse_dss. Indexed on bus (String), phase (integer),
then time (integer)

- TODO other reactive power specifications ?
- TODO handle vectors ("yearly" loads). Does the parse_dss handle redirects to txt files for loads?
- TODO LoadXfmrs.dss in 8500 Node ntwk: all loads are put on phases 1 or 2, 
    but actual MV phase could be 1, 2, or 3 -> have to map load phases through transformers?
"""
function dss_loads(d::Dict)
    P, Q = Dict{String, Dict{Int, Array{Real}}}(), Dict{String, Dict{Int, Array{Real}}}()
    for v in values(d["load"])
        bus = chop(v["bus1"], tail=length(v["bus1"])-findfirst('.', v["bus1"])+1)
        phases = collect(parse(Int,ph) for ph in split(v["bus1"][findfirst('.', v["bus1"])+1:end], "."))
        if !(bus in keys(P))
            P[bus] = Dict{Int, Array{Real}}()
            Q[bus] = Dict{Int, Array{Real}}()
        end
        if v["phases"] == 1 && get(v, "conn", "") != DELTA  # DELTA is a PMD Enum
            P[bus][phases[1]] = [v["kw"] * 1000]  
            if "kvar" in keys(v)
                Q[bus][phases[1]] = [v["kvar"] * 1000]
            elseif "pf" in keys(v)
                p = P[bus][phases[1]]
                Q[bus][phases[1]] = sqrt( (p/v["pf"])^2 - p^2 )
            else
                Q[bus][phases[1]] = 0.0
            end
        else  # split the load evenly across phases
            p = v["kw"] / length(phases) * 1000
            
            if "kvar" in keys(v)
                q = v["kvar"] / length(phases) * 1000
            elseif "pf" in keys(v)
                q = sqrt( (p/v["pf"])^2 - p^2 )
            else
                q = 0.0
            end
            for phs in phases
                P[bus][phs] = [p]
                Q[bus][phs] = [q]
            end
        end
    end
    return P, Q
end
