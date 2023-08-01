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
