


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
    init_inputs!(mg::MetaDiGraph; init_vs::Dict = Dict())

Use the `:load_sum_order` in `mg` to `init_inputs!` in the correct order, i.e. set the loads
at the leaf - substation connections as sums of all the loads (and the voltages at substations)
"""
function init_inputs!(mg::MetaDiGraph; init_vs::Dict = Dict())
    lso = get_prop(mg, :load_sum_order)
    ps = [mg[v, :p] for v in lso]
    init_inputs!(ps; init_vs = init_vs)
end


function set_inputs!(mg::MetaDiGraph)
    for v in get_prop(mg, :load_sum_order) # ~breadth first search of vertices
        if !(has_prop(mg, v, :m))
            throw(@error "MetaDiGraph must have model properties :m for all vertices to use set_inputs!")
        end
        # if v has inneighbors use their voltages at connections
        p_below = get_prop(mg, v, :p)
        for v_above in inneighbors(mg, v)
            m_above = get_prop(mg, v_above, :m)
            p_below.v0 = sqrt.(value.(m_above[:vsqrd][p_below.substation_bus, :])).data  # vector of time
        end
        # if v has outneighbors then use v_below's m[:Pj][p_below.substation_bus,:] * p_above.Sbase as v's Pload at the same bus
        p_above = get_prop(mg, v, :p)
        for v_below in outneighbors(mg, v)
            m_below = get_prop(mg, v_below, :m)
            p_below = get_prop(mg, v_below, :p)
            p_above.Pload[p_below.substation_bus] = value.(m_below[:Pj][p_below.substation_bus, :]).data * p_above.Sbase
            p_above.Qload[p_below.substation_bus] = value.(m_below[:Qj][p_below.substation_bus, :]).data * p_above.Sbase
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
    unique!(at_busses)
    mg = MetaDiGraph()
    # initial split
    p_above, p_below = BranchFlowModel.split_inputs(p, at_busses[1]);
    add_vertex!(mg, :p, p_above)
    add_vertex!(mg, :p, p_below)
    add_edge!(mg, 1, 2)
    set_indexing_prop!(mg, :p)

    for (i,b) in enumerate(at_busses[2:end])
        # find the vertex to split
        vertex = 0
        for v in vertices(mg)
            if b in mg[v, :p].busses
                vertex = v
                break
            end
        end
        p_above, p_below = BranchFlowModel.split_inputs(mg[vertex, :p], b);
        set_prop!(mg, vertex, :p, p_above)  # replace the already set :p, which preserves inneighbors
        add_vertex!(mg, :p, p_below)  # vertex i+2
        if !isempty(outdegree(mg, vertex))
            # already have edge(s) for vertex -> outneighbors(mg, vertex)
            # but now i+2 could be the parent for some of the outneighbors(mg, vertex)
            outns = copy(outneighbors(mg, vertex))
            for neighb in outns
                if !( mg[neighb, :p].substation_bus in mg[vertex, :p].busses )
                    # mv the edge to the new intermediate node
                    rem_edge!(mg, vertex, neighb)
                    add_edge!(mg, i+2, neighb)
                end
            end
        end
        add_edge!(mg, vertex, i+2)  # p_above -> p_below
    end
    # create the load_sum_order, a breadth first search from the leafs

    function recur_inneighbors(mg, vs, ins)
        for v in vs
            union!(ins, inneighbors(mg, v))
        end
        # do not combine these for loops (to keep order of `ins`)
        for v in vs
            invs = inneighbors(mg, v)
            return recur_inneighbors(mg, invs, ins)
        end
        return ins
    end
    leafvs = leaf_vertices(mg)
    set_prop!(mg, :load_sum_order, recur_inneighbors(mg, leafvs, leafvs))
    init_inputs!(mg)

    return mg
end



function get_diffs(mg::MetaDiGraph)
    pdiffs, qdiffs, vdiffs = Float64[], Float64[], Float64[]
    for v in get_prop(mg, :load_sum_order) # ~breadth first search of vertices
        # if v has inneighbors use their voltages at connections
        p_below = get_prop(mg, v, :p)
        for v_above in inneighbors(mg, v)
            m_above = get_prop(mg, v_above, :m)
            v_above = sqrt.(value.(m_above[:vsqrd][p_below.substation_bus, :]))  # vector of time
            push!(vdiffs, sum(abs.(p_below.v0 .- v_above)) / p_below.Ntimesteps)
        end

        m_above = get_prop(mg, v, :m)
        for v_below in outneighbors(mg, v)
            m_below = get_prop(mg, v_below, :m)
            b = get_prop(mg, v_below, :p).substation_bus
            push!(pdiffs, 
                sum(abs.(value.(m_below[:Pj][b,:]).data + value.(m_above[:Pj][b,:]).data)) / p_below.Ntimesteps
            )
            push!(qdiffs, 
                sum(abs.(value.(m_below[:Qj][b,:]).data + value.(m_above[:Qj][b,:]).data)) / p_below.Ntimesteps
            )
        end
    end
    return pdiffs, qdiffs, vdiffs
end