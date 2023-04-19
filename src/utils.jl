"""
    rij(i::AbstractString, j::AbstractString, p::Inputs{SinglePhase})

The per-unit resistance of line i->j
"""
function rij(i::AbstractString, j::AbstractString, p::Inputs{SinglePhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    rmatrix = p.Zdict[linecode]["rmatrix"] * linelength / p.Zbase
    return rmatrix[1]
end


"""
    xij(i::AbstractString, j::AbstractString, p::Inputs{SinglePhase})

The per-unit reacttance of line i->j
"""
function xij(i::AbstractString, j::AbstractString, p::Inputs{SinglePhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    xmatrix = p.Zdict[linecode]["xmatrix"] * linelength / p.Zbase
    return xmatrix[1]
end


function rij(i::AbstractString, j::AbstractString, p::Inputs{MultiPhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    rmatrix = p.Zdict[linecode]["rmatrix"] * linelength / p.Zbase
    return rmatrix
end


function xij(i::AbstractString, j::AbstractString, p::Inputs{MultiPhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    xmatrix = p.Zdict[linecode]["xmatrix"] * linelength / p.Zbase
    return xmatrix
end


# TODO test this function (esp. for matching openDSS results for IEEE13 system)
function zij(i::AbstractString, j::AbstractString, p::Inputs{MultiPhase})
    r = rij(i, j, p)
    x = xij(i, j, p)
    # need to fill out to 3x3 with zeros
    phases = sort(p.phases_into_bus[j])
    z = convert(Matrix{ComplexF64}, zeros(3,3))
    for (ii, phs1) in enumerate(phases), (jj, phs2) in enumerate(phases)
        z[phs1, phs2] = r[ii, jj] + x[ii, jj]*im
    end
    return z
end


"""
    remove_bus!(j::String, p::Inputs{SinglePhase})

Remove bus `j` in the line i->j->k from the model by making an equivalent line from busses i->k
"""
function remove_bus!(j::String, p::Inputs{SinglePhase})
    # get all the old values
    i, k = i_to_j(j, p)[1], j_to_k(j, p)[1]
    ij_idx, jk_idx = get_ij_idx(i, j, p), get_ij_idx(j, k, p)
    ij_len, jk_len = p.linelengths[ij_idx], p.linelengths[jk_idx]
    ij_linecode, jk_linecode = get_ijlinecode(i,j,p), get_ijlinecode(j,k,p)
    r_ij, x_ij, r_jk, x_jk = rij(i,j,p)*p.Zbase, xij(i,j,p)*p.Zbase, rij(j,k,p)*p.Zbase, xij(j,k,p)*p.Zbase
    # make the new values
    r_ik = r_ij + r_jk
    x_ik = x_ij + x_jk
    ik_len = ij_len + jk_len
    ik_linecode = ik_key = i * "-" * k
    ik_amps = minimum([p.Isquared_up_bounds[ij_linecode], p.Isquared_up_bounds[jk_linecode]])
    # delete the old values
    delete_edge_ij!(i, j, p)
    delete_edge_ij!(j, k, p)
    delete_bus_j!(j, p)
    # add the new values
    push!(p.edges, (i, k))
    push!(p.linecodes, ik_linecode)
    push!(p.phases, [1])
    push!(p.linelengths, ik_len)
    push!(p.edge_keys, ik_key)
    p.Zdict[ik_linecode] = Dict(
        "nphases" => 1,
        "name" => ik_linecode,
        "rmatrix" => [r_ik / ik_len],
        "xmatrix" => [x_ik / ik_len],
    )
    p.Isquared_up_bounds[ik_linecode] = ik_amps
end


"""
    remove_bus!(j::String, p::Inputs{MultiPhase})

Remove bus `j` in the line i->j->k from the model by making an equivalent line from busses i->k
"""
function remove_bus!(j::String, p::Inputs{MultiPhase})
    # get all the old values
    i, k = i_to_j(j, p)[1], j_to_k(j, p)[1]
    ij_idx, jk_idx = get_ij_idx(i, j, p), get_ij_idx(j, k, p)
    ij_len, jk_len = p.linelengths[ij_idx], p.linelengths[jk_idx]
    ij_linecode, jk_linecode = get_ijlinecode(i,j,p), get_ijlinecode(j,k,p)
    r_ij, x_ij, r_jk, x_jk = rij(i,j,p)*p.Zbase, xij(i,j,p)*p.Zbase, rij(j,k,p)*p.Zbase, xij(j,k,p)*p.Zbase
    phases = p.phases[ij_idx]
    # make the new values
    r_ik = r_ij .+ r_jk
    x_ik = x_ij .+ x_jk
    ik_len = ij_len + jk_len
    ik_linecode = ik_key = i * "-" * k
    ik_amps = minimum([p.Isquared_up_bounds[ij_linecode], p.Isquared_up_bounds[jk_linecode]])
    # delete the old values
    delete_edge_ij!(i, j, p)
    delete_edge_ij!(j, k, p)
    delete_bus_j!(j, p)
    # add the new values
    push!(p.edges, (i, k))
    push!(p.linecodes, ik_linecode)
    push!(p.phases, phases)
    push!(p.linelengths, ik_len)
    push!(p.edge_keys, ik_key)
    p.Zdict[ik_linecode] = Dict(
        "nphases" => length(phases),
        "name" => ik_linecode,
        "rmatrix" => r_ik ./ ik_len,
        "xmatrix" => x_ik ./ ik_len,
    )
    p.Isquared_up_bounds[ik_linecode] = ik_amps
end


"""
    combine_parallel_lines!(p::Inputs)

Combine any parallel single phase lines without loads on intermediate busses into one multiphase line

TODO this can probably be used for LinDistFlow too but then have to move rij and other
    utils.jl stuff to CommonOPF.jl (and distinguish between multiphase LinDistFlow and 
    BranchFlowModel Inputs)
"""
function combine_parallel_lines!(p::Inputs)
    g = make_graph(p.busses, p.edges)
    end_bs = busses_with_multiple_inneighbors(g)

    for b2 in end_bs
        ins = inneighbors(g, b2)
        start_bs = unique(
            next_bus_above_with_outdegree_more_than_one.(repeat([g], length(ins)), ins)
        )
        if length(start_bs) == 1 && typeof(start_bs[1]) == String
            # we have a start bus and end bus to merge lines (if none of the intermediate busses have loads)
            b1 = start_bs[1]
            paths = paths_between(g, b1, b2)
            check_paths(paths, p)
            # remove all the intermdiate busses s.t. we have two // lines from b1 to b2
            for path in paths
                for b in path
                    remove_bus!(b, p)
                end
            end
            # now we combine the two // lines into one
            ij_idxs = findall(t->(t[1]==b1 && t[2]==b2), p.edges)
            @assert length(ij_idxs) in [2,3] "Found more than three parallel lines between busses $b1 and $b2!"
            if length(ij_idxs) == 2
                # old values
                i1, i2 = ij_idxs
                len1, len2 = p.linelengths[i1], p.linelengths[i2]
                phases1, phases2 = p.phases[i1], p.phases[i2]
                linecode1, linecode2 = p.linecodes[i1], p.linecodes[i2]
                rmatrix1, rmatrix2 = p.Zdict[linecode1]["rmatrix"] .* len1, p.Zdict[linecode2]["rmatrix"] .* len2
                xmatrix1, xmatrix2 = p.Zdict[linecode1]["xmatrix"] .* len1, p.Zdict[linecode2]["xmatrix"] .* len2
                amps1, amps2 = p.Isquared_up_bounds[linecode1], p.Isquared_up_bounds[linecode2]
                # new values
                new_len = (len1 + len2) / 2
                new_linecode = "combined__" * linecode1
                new_rmatrix = Diagonal([rmatrix1[1], rmatrix2[1]]) ./ new_len
                new_xmatrix = Diagonal([xmatrix1[1], xmatrix2[1]]) ./ new_len
                # TODO add off diagonal r/x values from up/downstream lines?
                new_phases = [phases1[1], phases2[1]]

                # delete the old values
                delete_edge_index!(i1, p)
                ij_idxs = findall(t->(t[1]==b1 && t[2]==b2), p.edges)
                delete_edge_index!(ij_idxs[1], p)
                
                # add the new values
                push!(p.edges, (b1, b2))
                push!(p.linecodes, new_linecode)
                push!(p.phases, new_phases)
                push!(p.linelengths, new_len)
                push!(p.edge_keys, b1 * "-" * b2)
                p.phases_into_bus[b2] = new_phases
                p.Zdict[new_linecode] = Dict(
                    "nphases" => 2,
                    "name" => new_linecode,
                    "rmatrix" => new_rmatrix,
                    "xmatrix" => new_xmatrix,
                )
                p.Isquared_up_bounds[new_linecode] = (amps1 + amps2) / 2
            else  # 3 lines to combine
                # old values
                i1, i2, i3 = ij_idxs
                len1, len2, len3 = p.linelengths[i1], p.linelengths[i2], p.linelengths[i3]
                phases1, phases2, phases3 = p.phases[i1], p.phases[i2], p.phases[i3]
                linecode1, linecode2, linecode3 = p.linecodes[i1], p.linecodes[i2], p.linecodes[i3]
                rmatrix1, rmatrix2, rmatrix3 = p.Zdict[linecode1]["rmatrix"] .* len1, p.Zdict[linecode2]["rmatrix"] .* len2, p.Zdict[linecode3]["rmatrix"] .* len3
                xmatrix1, xmatrix2, xmatrix3 = p.Zdict[linecode1]["xmatrix"] .* len1, p.Zdict[linecode2]["xmatrix"] .* len2, p.Zdict[linecode3]["xmatrix"] .* len3
                amps1, amps2, amps3 = p.Isquared_up_bounds[linecode1], p.Isquared_up_bounds[linecode2], p.Isquared_up_bounds[linecode3]
                # new values
                new_len = (len1 + len2 + len3) / 3
                new_linecode = "combined__" * linecode1 *"__"* linecode2 *"__"* linecode3
                new_rmatrix = Diagonal([rmatrix1[1], rmatrix2[1], rmatrix3[1]]) ./ new_len
                new_xmatrix = Diagonal([xmatrix1[1], xmatrix2[1], xmatrix3[1]]) ./ new_len
                # TODO add off diagonal r/x values from up/downstream lines?
                new_phases = [phases1[1], phases2[1], phases3[1]]

                # delete the old values
                delete_edge_index!(i1, p)
                ij_idxs = findall(t->(t[1]==b1 && t[2]==b2), p.edges)
                delete_edge_index!(ij_idxs[1], p)
                ij_idxs = findall(t->(t[1]==b1 && t[2]==b2), p.edges)
                delete_edge_index!(ij_idxs[1], p)
                
                # add the new values
                push!(p.edges, (b1, b2))
                push!(p.linecodes, new_linecode)
                push!(p.phases, new_phases)
                push!(p.linelengths, new_len)
                push!(p.edge_keys, b1 * "-" * b2)
                p.phases_into_bus[b2] = new_phases
                p.Zdict[new_linecode] = Dict(
                    "nphases" => 3,
                    "name" => new_linecode,
                    "rmatrix" => new_rmatrix,
                    "xmatrix" => new_xmatrix,
                )
                p.Isquared_up_bounds[new_linecode] = (amps1 + amps2 + amps3) / 3
            end
            @info "Made new combined line between busses $b1 and $b2"
        end
    end
end



function get_constraints_by_variable_name(m::JuMP.AbstractModel, v::AbstractString)
    ac = ConstraintRef[]
    for tup in list_of_constraint_types(m)
        append!(ac, all_constraints(m, tup[1], tup[2]))
    end
    filter( cr -> occursin(v, string(cr)), ac )
end


"""
    get_load_bal_shadow_prices(m::JuMP.AbstractModel, p::Inputs)

create and return a dict indexed by bus and time for shadow prices
    (just real prices for now)
"""
function get_load_bal_shadow_prices(m::JuMP.AbstractModel, p::Inputs)
    @assert has_duals(m)
    d = Dict{String, Dict}(j => Dict{Int64, Float64}() for j in p.busses)

    for j in p.busses, t in 1:p.Ntimesteps
        d[j][t] = JuMP.shadow_price(m[:loadbalcons][j]["p"][t])
    end
    return d
end


"""
    cj(A)

short cut for conj(transpose(A))
"""
function cj(A)
    conj(transpose(A))
end


"""
    phi_ij(j::String, p::Inputs, M::AbstractMatrix)

Down-select the matrix M by the phase from i -> j
"""
function phi_ij(j::String, p::Inputs, M::AbstractMatrix)
    N = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
    for x in p.phases_into_bus[j], y in p.phases_into_bus[j]
        N[x,y] = M[x,y]
    end
    return N
end


"""
    check_unique_solution_conditions(p::Inputs)


report the maximum per-unit immpedance and load values.
See Chiang and Baran 2013:
A load flow solution with feasible voltage magnitude always exists and is unique when
1. V0 ≈ 1
2. loss values < 1
3. rpu, xpu << 1
"""
function check_unique_solution_conditions(p::Inputs)
    Rmax, Xmax = info_max_rpu_xpu(p)
    maxP, maxQ = info_max_Ppu_Qpu(p)

    if p.v0 > 1.09 || p.v0 < 0.91
        @warn "The substation voltage should be set close to one."
    end

    if maxP > 0.1 || maxQ > 0.1 
        @warn "The maximum load is greater than 0.1 * Sbase. You should probably increase Sbase."
    end
    if Rmax > 1e-2 || Xmax > 1e-2
        @warn "\nThe per unit impedance values should be much less than one but the max is >0.01 .\
            \nYou should increase Zbase by increasing Vbase."
    end
    nothing
end


"""
    reg_busses(p::Inputs)

All of the regulated busses, i.e. the second bus in the regulated edges
"""
function reg_busses(p::Inputs)
    getindex.(keys(p.regulators), 2)
end


function turn_ratio(p::Inputs, b::String)
    if !(b in reg_busses(p))
        throw(@error "Bus $b is not a regulated bus")
    end
    for (edge_tuple, d) in p.regulators
        if edge_tuple[2] == b
            return d[:turn_ratio]
        end
    end
end


function has_vreg(p::Inputs, b::String)
    for (edge_tuple, d) in p.regulators
        if edge_tuple[2] == b  && :vreg in keys(d)
            return true
        end
    end
    return false
end


function vreg(p::Inputs, b::String)
    if !(b in reg_busses(p))
        throw(@error "Bus $b is not a regulated bus")
    end
    for (edge_tuple, d) in p.regulators
        if edge_tuple[2] == b  && :vreg in keys(d)
            return d[:vreg]
        end
    end
    false
end


"""
    info_max_rpu_xpu(p::Inputs)

report the maximum per-unit resistance and reactance values of the lines.
It is important that the rpu ans xpu values be << 1. See Chiang and Baran 2013:
A load flow solution with feasible voltage magnitude always exists and is unique when
1. V0 ≈ 1
2. loss values < 1
3. rpu, xpu << 1
"""
function info_max_rpu_xpu(p::Inputs)
    Rmax = maximum([rij(i,j,p) for (i,j) in p.edges])
    Xmax = maximum([xij(i,j,p) for (i,j) in p.edges])
    @info("Max. Rpu: $Rmax   Max Xpu: $Xmax")
    return Rmax, Xmax
end
