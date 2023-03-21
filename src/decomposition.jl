


function leaf_busses(p::Inputs)
    leafs = String[]
    for j in p.busses
        if !isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))
            push!(leafs, j)
        end
    end
    return leafs
end


function leaf_vertices(mg::MetaDiGraph)
    leafs = Int64[]
    for v in vertices(mg)
        if outdegree(mg, v) == 0
            push!(leafs, v)
        end
    end
    return leafs
end

"""
    init_inputs!(ps::Vector{Inputs{BranchFlowModel.SinglePhase}}; init_vs::Dict = Dict())

Set the load on the upstream leaf noades equal to the sum of all the loads in the
downstream inputs. It is important that the order of `ps` is from leaf branches to trunk branches
so that the sums of loads take into account all downstream sub-trees.

if `init_vs` is provided, the `p.v0` is set for the `Input` with its `p.substation_bus` equal 
    to the key in `init_vs`
```julia
init_vs = Dict(
    "sub_bus_1" => 0.98
)

for p in ps
    if p.substation_bus in keys(init_vs)
        p.v0 = init_vs[p.substation_bus]
    end
end
```
"""
function init_inputs!(ps::Vector{Inputs{BranchFlowModel.SinglePhase}}; init_vs::Dict = Dict())
    for p in ps
        if p.substation_bus in keys(init_vs)
            p.v0 = init_vs[p.substation_bus]
        end
        leafs = leaf_busses(p)
        for pp in ps
            if pp.substation_bus in leafs
                p.Pload[pp.substation_bus] = sum( values(pp.Pload) )
                p.Qload[pp.substation_bus] = sum( values(pp.Qload) )
            end
        end
    end
    true
end


"""
    split_at_busses(p::Inputs{BranchFlowModel.SinglePhase}, at_busses::Vector{String})

Split `Inputs` using the `at_busses`

TODO? maybe add a `mode` for how to split, and a `max_vars` parameter

returns MetaDiGraph with vertex properties `:p` containing `Input` for the sub-graphs.
For example `mg[2, :p]` is the `Input` at the second vertex of the graph created by splitting 
the network via the `at_busses`.
"""
function split_at_busses(p::Inputs{BranchFlowModel.SinglePhase}, at_busses::Vector{String})
    mg = MetaDiGraph()
    g = make_graph(p.busses, p.edges)
    # initial split
    p_above, p_below = BranchFlowModel.split_inputs(p, at_busses[1], g);
    add_vertex!(mg, :p, p_above)
    add_vertex!(mg, :p, p_below)
    add_edge!(mg, 1, 2)
    set_indexing_prop!(mg, :p)

    for (i,b) in enumerate(at_busses[2:end])
        p_above, p_below = BranchFlowModel.split_inputs(p, b, g);
        # find the vertex that was just split
        vertex = 0
        for v in vertices(mg)
            if b in mg[v, :p].busses
                vertex = v
                break
            end
        end
        set_prop!(mg, vertex, :p, p_above)  # replace the already set :p
        add_vertex!(mg, :p, p_below)
        add_edge!(mg, vertex, i+2)
    end
    # create the load_sum_order, a breadth first search from the leafs

    function recur_inneighbors(mg, vs, ins)
        for v in vs
            union!(ins, inneighbors(mg, v))
        end
        for v in vs
            invs = inneighbors(mg, v)
            return recur_inneighbors(mg, invs, ins)
        end
        return ins
    end
    leafvs = leaf_vertices(mg)
    set_prop!(mg, :load_sum_order, recur_inneighbors(mg, leafvs, leafvs))
    return mg
end


"""
    init_inputs!(mg::MetaDiGraph; init_vs::Dict = Dict())

Use the `:load_sum_order` in `mg` to `init_inputs!` in the correct order, i.e. set the loads
at the leaf - substation connections as sums of all the loads
"""
function init_inputs!(mg::MetaDiGraph; init_vs::Dict = Dict())
    lso = get_prop(mg, :load_sum_order)
    ps = [mg[v, :p] for v in lso]
    init_inputs!(ps; init_vs = init_vs)
end