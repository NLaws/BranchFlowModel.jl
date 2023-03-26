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
    fill_transformer_vals!(d::Dict, p::Inputs)

Fill in necessary values for parsing transformers in the dict from parse_dss
"""
function fill_transformer_vals!(d::Dict, Sbase::Real, Vbase::Real)
    d["%r"]   = get(d, "%r",   0.01)
    d["%r_2"] = get(d, "%r_2", 0.01)

    if "kvs" in keys(d)
        d["kv"], d["kv_2"] = d["kvs"]
    else
        d["kv"]   = get(d, "kv", Vbase / 1000)
        d["kv_2"] = get(d, "kv_2", Vbase / 1000)
    end

    if "buses" in keys(d)
        d["bus"], d["bus_2"] = d["buses"]
    else
        d["bus"]   = get(d, "bus", nothing)
        d["bus_2"] = get(d, "bus_2", nothing)
    end

    if "kvas" in keys(d)
        d["kva"], d["kva_2"] = d["kvas"]
    else
        d["kva"]   = get(d, "kva", Sbase / 1000)
        d["kva_2"] = get(d, "kva_2", Sbase / 1000)
    end

    d["xlt"] = get(d, "xlt", 1)
    d["xht"] = get(d, "xlt", 1)
    d["xhl"] = get(d, "xlt", 1)
end


"""
    dss_dict_to_arrays(d::Dict)

Parse the dict from PowerModelsDistribution.parse_dss into values needed for LinDistFlow

TODO need a standardize format for the output of parse_dss rather than following all the ways
    of defining things in openDSS (for example busses can be defined as an array or individual 
    call outs for some objects).
"""
function dss_dict_to_arrays(d::Dict, Sbase::Real, Vbase::Real)
    # TODO allocate empty arrays with number of lines
    # TODO separate this method into sub-methods, generally parse components separately
    edges = Tuple[]
    phases = Vector[]
    linecodes = String[]
    linelengths = Float64[]
    Isquared_up_bounds = Dict{String, Float64}()

    if !("linecode" in keys(d))
        d["linecode"] = Dict{String, Any}()  # we add it b/c we use it for storing R/X values
    end

    # some reuseable stuff
    function get_b1_b2_phs(v::Dict)
        if occursin(".", v["bus1"])  # have to account for .1.2 phases for example
            b1 = chop(v["bus1"], tail=length(v["bus1"])-findfirst('.', v["bus1"])+1)
            phs = sort!(collect(parse(Int,ph) for ph in split(v["bus1"][findfirst('.', v["bus1"])+1:end], ".")))
        else  # default to 3 phases
            b1 = v["bus1"]
            if v["phases"] == 1
                phs = [1]
            else
                phs = [1,2,3]
            end
        end

        if occursin(".", v["bus2"])
            b2 = chop(v["bus2"], tail=length(v["bus2"])-findfirst('.', v["bus2"])+1)
        else
            b2 = v["bus2"]
        end
        return b1, b2, phs
    end

    switches_to_check = String[]  # only add switches if they are not in parallel with a line
    for (k,v) in d["line"]  # Line dict includes switches
        if "switch" in keys(v) && v["switch"] == true
            push!(switches_to_check, k)
            continue
        end
        try
            b1, b2, phs = get_b1_b2_phs(v)
            if b2 in tails(edges)
                @warn "Not adding line $k ($b1 to $b2) because there is already an edge into $b2"
                continue
            end
            push!(edges, (b1, b2))
            push!(phases, phs)

            if "linecode" in keys(v)
                push!(linecodes, v["linecode"])
            else  # assume positive sequence
                v["linecode"] = v["name"]
                push!(linecodes, v["linecode"])
                d["linecode"][v["name"]] = Dict{String, Any}(
                    "rmatrix" => [v["r1"]],
                    "xmatrix" => [v["x1"]],
                    "nphases" => 1,
                    "name" => v["name"]
                )
            end

            # TODO ratings could be in linecode dict too
            # TODO assuming that there are linecodes, should converge on consistent keys for lines
            Isquared_up_bounds[v["linecode"]] = DEFAULT_AMP_LIMIT^2
            if "normamps" in keys(v) && !(v["normamps"] ≈ 0)  # assuming lowercase keys
                Isquared_up_bounds[v["linecode"]] = v["normamps"]^2
            elseif "emergamps" in keys(v) && !(v["emergamps"] ≈ 0)
                Isquared_up_bounds[v["linecode"]] = v["emergamps"]^2
            end

            # TODO handle scaling of lengths and R/X values
            # for now just make sure the linecode and line values are in consistent units
            # and BEWARE PowerModelsDistribution will scale values from openDSS!
            push!(linelengths, v["length"]) 
        catch
            @warn("Unable to parse line $(k) when processing OpenDSS model.")
        end
    end

    for k in switches_to_check
        v = d["line"][k]
        try
            # need to connect busses over switch
            b1, b2, phs = get_b1_b2_phs(v)
            if (b1,b2) in edges
                @warn "Not adding switch $k because there is already an edge from $b1 to $b2)"
                continue
            end
            linecode = "switch" * v["name"]
            push!(edges, (b1, b2))
            push!(linecodes, linecode)
            push!(linelengths, get(v, "length", 1.0))
            push!(phases, phs)

            Isquared_up_bounds[linecode] = DEFAULT_AMP_LIMIT^2
            if "normamps" in keys(v) && !(v["normamps"] ≈ 0)  # assuming lowercase keys
                Isquared_up_bounds[linecode] = v["normamps"]^2
            elseif "emergamps" in keys(v) && !(v["emergamps"] ≈ 0)
                Isquared_up_bounds[linecode] = v["emergamps"]^2
            end

            d["linecode"][linecode] = Dict(
                "nphases" => length(phs),
                "rmatrix" => Diagonal(ones(3)) * get(v, "r1", 0.00001),
                "xmatrix" => Diagonal(ones(3)) *  get(v, "x1", 0.00001),
            )
        catch e
            @warn("Unable to parse switch $(k) when processing OpenDSS model.")
            println(e)
        end
    end

    # make phases_into_bus to infer transformer phases
    phases_into_bus = Dict(k=>v for (k,v) in zip(tails(edges), phases))

    # TODO it is possible to have a transformer -> transformer, which could lead to failed
    # parsing
    transformers_to_try_again = String[]
    trfxs_with_regs = [innerd["transformer"] for (k, innerd) in get(d, "regcontrol", Dict())]
    regulators = Dict()
    for (k,v) in get(d, "transformer", Dict())
        try
            # need to connect busses over transformers
            fill_transformer_vals!(v, Sbase, Vbase)  # need Sbase, Vbase
            b1, b2 = v["bus"], v["bus_2"]

            if !(b1 in tails(edges)) && !(b2 in heads(edges))
                # this transformer does not connect anything so we ignore it
                @warn("Not parsing transformer $k between $b1 and $b2
                    because it does not have an edge in to or out of it.")
                continue
            end

            if b1 in tails(edges)
                phs = phases_into_bus[b1]
            else
                push!(transformers_to_try_again, k)
                continue
            end

            nwindings = v["windings"]
            if nwindings != 2
                @warn("Parsing a $nwindings winding transformer as a 2 winding transformer.")
            end

            R1 = v["%r"] / 100 * v["kv"]^2 * 1000 / v["kva"]
            R2 = v["%r_2"] / 100 * v["kv"]^2 * 1000 / v["kva"]
            R = R1 + R2
            X = (v["xhl"] + v["xlt"] + v["xht"]) / 100 * v["kv"]^2 * 1000 / v["kva"]
            # openDSS Manual says "Always use the kVA base of the first winding for entering impedances. 
            # Impedance values are entered in percent."

            linecode = v["name"]
            push!(edges, (b1, b2))
            push!(linecodes, linecode)
            push!(linelengths, 1.0)
            push!(phases, phs)
            phases_into_bus[b2] = phs

            Isquared_up_bounds[linecode] = DEFAULT_AMP_LIMIT^2

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
            # TODO parse turn ratios
            if linecode in trfxs_with_regs
                regulators[(b1,b2)] = Dict(:turn_ratio => 1.0)
            end

        catch e
            @warn("Unable to parse transformer $(k) when processing OpenDSS model.")
            println(e)
        end
    end

    for k in transformers_to_try_again
        v = d["transformer"][k]
        b1, b2 = v["bus"], v["bus_2"]

        if b1 in tails(edges)
            phs = phases_into_bus[b1]
        else  # should be feeder head transformer but TODO add checks
            if all(length(phz) == 1 for phz in phases)
                phs = [1]
            else
                phs = [1,2,3]
            end
        end
        nwindings = get(v, "windings", 2)
        if nwindings != 2
            @warn("Parsing a $nwindings winding transformer as a 2 winding transformer.")
        end

        R1 = v["%r"] / 100 * v["kv"]^2 * 1000 / v["kva"]
        R2 = v["%r_2"] / 100 * v["kv"]^2 * 1000 / v["kva"]
        R = R1 + R2
        X = (v["xhl"] + v["xlt"] + v["xht"]) / 100 * v["kv"]^2 * 1000 / v["kva"]
        # openDSS Manual says "Always use the kVA base of the first winding for entering impedances. 
        # Impedance values are entered in percent."

        linecode = v["name"]
        push!(edges, (b1, b2))
        push!(linecodes, linecode)
        push!(linelengths, 1.0)
        push!(phases, phs)
        phases_into_bus[b2] = phs

        Isquared_up_bounds[linecode] = DEFAULT_AMP_LIMIT^2

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
        # TODO parse turn ratios
        if linecode in trfxs_with_regs
            regulators[(b1,b2)] = Dict(:turn_ratio => 1.0)
        end
    end

    return edges, linecodes, linelengths, d["linecode"], phases, Isquared_up_bounds, regulators
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

        if occursin(".", v["bus1"])
            bus = chop(v["bus1"], tail=length(v["bus1"])-findfirst('.', v["bus1"])+1)
            phases = collect(parse(Int,ph) for ph in split(v["bus1"][findfirst('.', v["bus1"])+1:end], "."))
        else
            bus = v["bus1"]
            phases = [1]
        end
        if !(bus in keys(P))
            P[bus] = Dict{Int, Array{Real}}(phases[1] => [0.0])
            Q[bus] = Dict{Int, Array{Real}}(phases[1] => [0.0])
        end
        if v["phases"] == 1 && get(v, "conn", "") != DELTA  # DELTA is a PMD Enum
            phs = phases[1]
            if phs in keys(P[bus])
                P[bus][phs][1] += v["kw"] * 1000
            else
                P[bus][phs] = [v["kw"] * 1000]
            end
            if "kvar" in keys(v)
                if phs in keys(Q[bus])
                    Q[bus][phs][1] += v["kvar"] * 1000
                else
                    Q[bus][phs] = [v["kvar"] * 1000]
                end
            elseif "pf" in keys(v)
                p = P[bus][phs][1]
                if phs in keys(Q[bus])
                    Q[bus][phs][1] += sqrt( (p/v["pf"])^2 - p^2 )
                else
                    Q[bus][phs] = [sqrt( (p/v["pf"])^2 - p^2 )]
                end
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
                if phs in keys(P[bus])
                    P[bus][phs][1] += p
                else
                    P[bus][phs] = [p]
                end
                if phs in keys(Q[bus])
                    Q[bus][phs][1] += q
                else
                    Q[bus][phs] = [q]
                end
            end
        end
    end
    return P, Q
end
