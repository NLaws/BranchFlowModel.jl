"""
    delete_edge_ij!(i::String, j::String, p::Inputs{SinglePhase})

delete edge `(i, j)` from
- p.edges
- p.phases
- p.linelengths
- p.edge_keys
- p.Isquared_up_bounds

NOTE do not delete!(p.Zdict, ij_linecode) nor delete!(p.Isquared_up_bounds, ij_linecode) 
because anything indexed on linecodes can be used for multiple lines
"""
function delete_edge_ij!(i::String, j::String, p::Inputs{SinglePhase})
    idx = get_ij_idx(i, j, p)
    deleteat!(p.edges,       idx)
    deleteat!(p.linecodes,   idx)
    deleteat!(p.phases,      idx)
    deleteat!(p.linelengths, idx)
    deleteat!(p.edge_keys,   idx)
    true
end


"""
    delete_bus_j!(j::String, p::Inputs{SinglePhase})

Remove bus `j` from `p.busses`
"""
function delete_bus_j!(j::String, p::Inputs{SinglePhase})
    p.busses = setdiff(p.busses, [j])
    if j in keys(p.Pload)
        delete!(p.Pload, j)
    end
    if j in keys(p.Qload)
        delete!(p.Qload, j)
    end
    true
end


"""
    reduce_tree!(p::Inputs{SinglePhase})

combine any line sets with intermediate busses that have indegree == outdegree == 1
and is not a load bus into a single line
"""
function reduce_tree!(p::Inputs{SinglePhase}; keep_regulator_busses=true)
    # TODO make graph once in Inputs ?
    g = make_graph(p.busses, p.edges)
    int_bus_map = get_prop(g, :int_bus_map)
    reducable_buses = String[]
    load_buses = Set(vcat(collect(keys(p.Pload)), collect(keys(p.Qload))))
    if keep_regulator_busses
        for bs in keys(p.regulators)  # tuple keys of bus pairs, i.e. edges
            push!(load_buses, bs...)
        end
    end
    for v in vertices(g)
        if indegree(g, v) == outdegree(g, v) == 1 && !(int_bus_map[v] in load_buses)
            push!(reducable_buses, int_bus_map[v])
        end
    end
    @debug("Removing the following busses: \n$reducable_buses")
    # replace two lines with one
    for j in reducable_buses
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
    @info("Removed $(length(reducable_buses)) busses.")
end


function trim_tree_once!(p::Inputs{SinglePhase})
    trimmable_busses = setdiff(leaf_busses(p), union(keys(p.Pload), keys(p.Qload)))
    if isempty(trimmable_busses) return false end
    trimmable_edges = Tuple[]
    for j in trimmable_busses
        for i in i_to_j(j, p)
            push!(trimmable_edges, (i,j))
        end
    end
    @debug("Deleting the following edges from the Inputs:")
    for edge in trimmable_edges @debug(edge) end
    for (i,j) in trimmable_edges
        delete_edge_ij!(i, j, p)
        delete_bus_j!(j, p)
    end
    true
end


function trim_tree!(p::Inputs{SinglePhase})
    n_edges_before = length(p.edges)
    trimming = trim_tree_once!(p)
    while trimming
        trimming = trim_tree_once!(p)
    end
    n_edges_after = length(p.edges)
    @info("Removed $(n_edges_before - n_edges_after) edges.")
    true
end


function make_sub_inputs(
    p::Inputs{SinglePhase}, 
    edges_to_delete::Vector, 
    busses_to_delete::Vector{String}
    )
    pc = deepcopy(p)
    for e in edges_to_delete
        delete_edge_ij!(e[1], e[2], pc)
    end
    for b in busses_to_delete
        delete_bus_j!(b, pc)
    end
    return pc
end


"""
    split_inputs(p::Inputs{SinglePhase}, bus::String, g::SimpleDiGraph)

Split inputs into one graph for everything above `bus` and one graph for everything
    below `bus`.
"""
function split_inputs(p::Inputs{SinglePhase}, bus::String)
    g = make_graph(p.busses, p.edges)
    in_buses = collect(all_inneighbors(g, bus, String[]))
    out_buses = collect(all_outneighbors(g, bus, String[], String[]))
    # in/out_buses do not have bus, but sub_busses does have bus
    # we want to keep bus in both Inputs

    sub_busses, sub_edges = induced_subgraph(g, vcat(out_buses, bus))
    p_above = make_sub_inputs(p, sub_edges, out_buses)

    sub_busses, sub_edges = induced_subgraph(g, vcat(in_buses, bus))
    p_below = make_sub_inputs(p, sub_edges, in_buses)
    p_below.substation_bus = bus

    return p_above, p_below
end



"""
    split_inputs(p::Inputs{SinglePhase}, bus::String, out_buses::Vector{String})

Split `p` into `p_above` and `p_below` where `p_below` has only `out_buses` and `p_above` has
`union( [bus], setdiff(p.busses, out_buses) )`.

Note that `out_buses` must contain `bus`
"""
function split_inputs(p::Inputs{SinglePhase}, bus::String, out_buses::Vector{String})
    if !(bus in out_buses)
        throw(@error "Cannot split inputs: bus is not in out_buses.")
    end
    g = make_graph(p.busses, p.edges)
    in_buses = setdiff(p.busses, out_buses)
    # in_buses does not have bus, but sub_busses does have bus
    # we want to keep bus in both p_above and p_below

    sub_busses, edges_to_remove = induced_subgraph(g, out_buses)
    # NOTE edges_to_remove does not include edges going out of out_buses
    p_above = make_sub_inputs(p, edges_to_remove, setdiff(out_buses, [bus]))

    sub_busses, edges_to_remove = induced_subgraph(g, vcat(in_buses, bus))
    p_below = make_sub_inputs(p, edges_to_remove, in_buses)
    p_below.substation_bus = bus

    return p_above, p_below
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


function info_max_Ppu_Qpu(p::Inputs)
    maxP = maximum(maximum.(values(p.Pload))) / p.Sbase
    maxQ = maximum(maximum.(values(p.Qload))) / p.Sbase
    @info("Max. Ppu: $maxP   Max Qpu: $maxQ")
    return maxP, maxQ
end