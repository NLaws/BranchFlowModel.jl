"""
    leaf_busses(p::Inputs)

returns `Vector{String}` containing all of the leaf busses in `p.busses`
"""
function leaf_busses(p::Inputs)
    leafs = String[]
    for j in p.busses
        if !isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))
            push!(leafs, j)
        end
    end
    return leafs
end


"""
    connecting_busses(mg::MetaDiGraph, v)

return Dict with keys for out-vertices of `v` and values for their substation_bus
"""
function connecting_busses(mg::MetaDiGraph, v)
    cbusses = Dict()
    for neighb in  outneighbors(mg, v)
        if mg[neighb, :p].substation_bus in mg[v, :p].busses
            cbusses[neighb] = mg[neighb, :p].substation_bus
        end
    end
    return cbusses
end


"""
    leaf_vertices(mg::MetaDiGraph)

returns `Vector{Int64}` containing all of the leaf vertices in `mg`
"""
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
                if isempty(pp.Pload)
                    @warn "The Pload dictionary is empty for the sub network with head bus $(pp.substation_bus).\
                    \nThis indicates that there are probably lines with no loads on the ends in the larger network.\
                    \nTry using trim_tree! to eliminate the deadend branches."
                else
                    p.Pload[pp.substation_bus] = sum( values(pp.Pload) )
                end
                if isempty(pp.Qload)
                    @warn "The Qload dictionary is empty for the sub network with head bus $(pp.substation_bus).\
                    \nThis indicates that there are probably lines with no loads on the ends in the larger network.\
                    \nTry using trim_tree! to eliminate the deadend branches."
                else
                    p.Qload[pp.substation_bus] = sum( values(pp.Qload) )
                end
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


"""
    set_inputs!(mg::MetaDiGraph; α::Float64=0.0)

Set the shared values in each subgraph / vertex of mg:
1. set the current vertex's v0 to its inneighbor's voltage
2. set the current vertex P/Qload to the outneighbors' substation_bus loads
"""
function set_inputs!(mg::MetaDiGraph; α::R=0.0) where R <: Real
    for v in get_prop(mg, :load_sum_order) # ~breadth first search of vertices
        # if v has inneighbors use their voltages at connections
        p_below = get_prop(mg, v, :p)
        for v_above in inneighbors(mg, v)
            m_above = get_prop(mg, v_above, :m)
            if α == 0.0
                p_below.v0 = sqrt.(value.(m_above[:vsqrd][p_below.substation_bus, :])).data  # vector of time
            else  # use weighted average of new and old values
                m_below = get_prop(mg, v, :m)
                v_kp1 = sqrt.(value.(m_above[:vsqrd][p_below.substation_bus, :])).data
                v_k   = sqrt.(value.(m_below[:vsqrd][p_below.substation_bus, :])).data
                p_below.v0 = (v_kp1 + α .* v_k) ./ (1 + α)
            end
        end
        # if v has outneighbors then use v_below's m[:Pj][p_below.substation_bus,:] * p_above.Sbase as v's Pload at the same bus
        # have set the loads from deepest vertices to shallowest to get correct sums
        p_above = get_prop(mg, v, :p)
        for v_below in outneighbors(mg, v)
            m_below = get_prop(mg, v_below, :m)
            p_below = get_prop(mg, v_below, :p)
            if α == 0.0
                p_above.Pload[p_below.substation_bus] = value.(m_below[:Pj][p_below.substation_bus, :]).data * p_above.Sbase
                p_above.Qload[p_below.substation_bus] = value.(m_below[:Qj][p_below.substation_bus, :]).data * p_above.Sbase
            else
                p_kp1 = value.(m_below[:Pj][p_below.substation_bus, :]).data * p_above.Sbase
                q_kp1 = value.(m_below[:Qj][p_below.substation_bus, :]).data * p_above.Sbase
                p_k = copy(p_above.Pload[p_below.substation_bus])
                q_k = copy(p_above.Qload[p_below.substation_bus])
                p_above.Pload[p_below.substation_bus] = (p_kp1 + α .* p_k) ./ (1 + α)
                p_above.Qload[p_below.substation_bus] = (q_kp1 + α .* q_k) ./ (1 + α)
            end
        end
    end
    true
end


"""
    breadth_first_nodes_from_leafs(mg::AbstractGraph)

returns a `Vector{Int64}` for the vertices in `mg` from the leafs to the source.
"""
function breadth_first_nodes_from_leafs(mg::AbstractGraph)
    order = leaf_vertices(mg)  # initialize with all leafs
    nodes_left = Set(setdiff(vertices(mg), order))
    while !isempty(nodes_left)
        for n in order
            inn = inneighbors(mg, n)  # only one per node except for trunk
            if length(inn) == 1
                inn = inn[1]
                if all(outn in order for outn in outneighbors(mg, inn))
                    # then inn can join the order b/c it has no nodes below it in nodes_left
                    union!(order, [inn])
                    delete!(nodes_left, inn)
                elseif isempty(inn) && length(nodes_left) == 1
                    # trunk is last
                    union!(order, [inn])
                    delete!(nodes_left, inn)
                end
            end
        end
    end
    return order
end


"""
    split_at_busses(p::Inputs{BranchFlowModel.SinglePhase}, at_busses::Vector{String})

Split `Inputs` using the `at_busses`

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
    vs, depths = vertices_from_deepest_to_source(mg, 1)
    set_prop!(mg, :load_sum_order, vs)
    init_inputs!(mg)
    if mg.graph.ne != length(mg.vprops) - 1
        @warn "The MetaDiGraph created is not a tree."
    end

    return mg
end


"""
    split_at_busses(p::Inputs{BranchFlowModel.SinglePhase}, at_busses::Vector{String}, with_busses::Vector{Vector{String}})

Split up `p` using the `at_busses` as each new `substation_bus` and containing the corresponding `with_busses`.
The `at_busses` and `with_busses` can be determined using `splitting_busses`.

NOTE: this variation of splt_at_busses allows for more than two splits at the same bus; whereas the other
implementation of split_at_busses only splits the network into two parts for everything above and
everything below a splitting bus.
"""
function split_at_busses(p::Inputs{BranchFlowModel.SinglePhase}, at_busses::Vector{String}, with_busses::Vector{Vector}; add_connections=true)
    unique!(at_busses)
    mg = MetaDiGraph()
    if add_connections
        with_busses = connect_subgraphs_at_busses(p, at_busses, with_busses)
    end
    # initial split
    p_above, p_below = BranchFlowModel.split_inputs(p, at_busses[1], with_busses[1]);
    add_vertex!(mg, :p, p_above)
    add_vertex!(mg, :p, p_below)
    add_edge!(mg, 1, 2)
    set_indexing_prop!(mg, :p)

    for (i, (b, sub_bs)) in enumerate(zip(at_busses[2:end], with_busses[2:end]))
        # find the vertex to split
        vertex = 0
        for v in vertices(mg)
            if b in mg[v, :p].busses
                vertex = v
                break
            end
        end
        p_above, p_below = BranchFlowModel.split_inputs(mg[vertex, :p], b, sub_bs);
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
    vs, depths = vertices_from_deepest_to_source(mg, 1)
    set_prop!(mg, :load_sum_order, vs)
    init_inputs!(mg)
    if mg.graph.ne != length(mg.vprops) - 1
        @warn "The MetaDiGraph created is not a tree."
    end

    return mg
end



"""
    get_diffs(mg::MetaDiGraph)

Uses the JuMP Models stored in mg[:m] to calculate the difference between Pj, Qj, and |v| at every
leaf/substation connection. 

returns three `Float64[]`
"""
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


"""
    splitting_busses(p::Inputs{BranchFlowModel.SinglePhase}, source::String; threshold::Int64=10)

Determine the busses to split a tree graph on by searching upward from the deepest leafs first
and gathering the nearest busses until threshold is met for each subgraph.

Returns a `Vector{String}` for the bus names and `Vector{Vector{String}}` for the corresponding
busses within each sub-graph.

!!! note
    It is not enough to have only the splitting busses to obey the `max_busses` limit because
    one must also know which sub branches to take from each splitting bus. In other words, we also
    need all the busses within each subgraph to split properly. For example, if a splitting
    bus has two sub branches then obeying the max_busses limit can require only including one
    sub branch out of the splitting bus. To know which branch to take we can use the other busses
    in the sub graph.
"""
function splitting_busses(p::Inputs{BranchFlowModel.SinglePhase}, source::String; max_busses::Int64=10)
    g = make_graph(p.busses, p.edges)
    bs, depths = busses_from_deepest_to_source(g, source)
    splitting_bs = String[]  # head nodes of all the subgraphs
    subgraph_bs = Vector[]
    # iterate until bs is empty, taking out busses as they are added to subgraphs
    subg_bs = String[]
    bs_parsed = String[]
    while !isempty(bs)
        b = popfirst!(bs)
        push!(subg_bs, b)
        ins = [b] # first check for any out neighbors of b
        while length(ins) == 1  # moving up tree from b in this loop
            inb = ins[1]
            # outns includes every bus below inb, 
            # excluding any branches that start with a bus in bs_parsed
            outns = all_outneighbors(g, inb, String[], bs_parsed)
            setdiff!(outns, bs_parsed)  # just in case
            new_subg_bs = unique(vcat([inb], outns, subg_bs))

            if length(new_subg_bs) > max_busses || isempty(bs)
                # addition of busses would increase busses in subgraph beyond max_busses
                # so we split at b and start a new subgraph
                push!(splitting_bs, b)
                push!(bs_parsed, subg_bs...)
                push!(subgraph_bs, subg_bs)
                bs = setdiff(bs, subg_bs)
                subg_bs = String[]
                break  # inner loop
            end
            # else continue going up tree
            subg_bs = new_subg_bs
            bs = setdiff(bs, subg_bs)
            ins = inneighbors(g, inb)  # go up another level
            b = inb
        end
    end
    if source in splitting_bs
        return setdiff(splitting_bs, [source]), subgraph_bs[1:end-1]
    end
    return splitting_bs, subgraph_bs
    # the last bus in splitting_bs is the source, which is not really a splitting bus
end


"""
    connect_subgraphs_at_busses(p::Inputs{BranchFlowModel.SinglePhase}, at_busses::Vector{String}, subgraphs::Vector{Vector})

The splitting_busses algorithm does not include over laps in subgraphs.
But, we want overlaps at the splitting busses for solving the decomposed branch flow model.
So here we add the overlapping splitting busses to each sub graph.
"""
function connect_subgraphs_at_busses(p::Inputs{BranchFlowModel.SinglePhase}, at_busses::Vector{String}, subgraphs::Vector{Vector})
    g = make_graph(p.busses, p.edges)
    new_subgs = deepcopy(subgraphs)
    for (i, subgraph) in enumerate(subgraphs)
        for b in subgraph
            bs_to_add = intersect(outneighbors(g, b), at_busses)
            if !isempty(bs_to_add)
                for ba in bs_to_add
                    if !(ba in new_subgs[i])
                        push!(new_subgs[i], ba)
                    end
                end
            end
        end
    end
    return new_subgs
end


"""
    metagraph_to_json(mg::MetaDiGraph, filename::String)

Dump the vertices, edges, and the busses and lines for each vertex in the metagraph to JSON
"""
function metagraph_to_json(mg::MetaDiGraph, filename::String)
    vs, order = vertices_from_deepest_to_source(mg, 1)
    depths = Dict(k => v for (k,v) in zip(vs, order))
    d = Dict()
    d["nodes"] = collect(vertices(mg))
    d["edges"] = [(e.src, e.dst) for e in collect(edges(mg))]
    d["depths"] = depths
    for v in vertices(mg)
        d[v] = Dict()
        d[v]["busses"] = mg[v, :p].busses
        d[v]["lines"] = mg[v, :p].edges
    end
    open(filename * ".json", "w") do f
        JSON.print(f, d, 2)
    end
end


"""
    build_metagraph(p::Inputs{BranchFlowModel.SinglePhase}, source::String; max_busses::Int64=10)

return MetaDiGraph with `:p` property set for every vertex by splitting the `Inputs` via `splitting_busses`
"""
function build_metagraph(p::Inputs{BranchFlowModel.SinglePhase}, source::String; max_busses::Int64=10)
    splitting_bs, subgraphs = splitting_busses(p, source; max_busses=max_busses)  
    split_at_busses(p, splitting_bs, subgraphs)
end


"""
    solve_metagraph!(mg::MetaDiGraph, builder::Function, tol::T; α::T=0.5, verbose=false) where T <: Real

Given a MetaDiGraph and a JuMP Model `builder` method iteratively solve the models until the `tol` is 
met for the differences provided by `BranchFlowModel.get_diffs`. 

The `builder` must accept only one argument of type `BranchFlowModel.AbstractInputs` that returns 
a `JuMP.AbstractModel`. Each model returned from the `builder` is stored as an `:m` property in 
each vertex of `mg`.

!!! note 
    `tol` is compared to the maximum absolute value of all the p, q, and v differences.
"""
function solve_metagraph!(mg::MetaDiGraph, builder::Function, tol::R; α::T=0.5, verbose=false) where {T <: Real, R <: Real}
    init_inputs!(mg)
    diff = abs(tol) * 10
    i = 0
    while diff > abs(tol)
        for (p, vertex) in mg[:p]
            m = builder(p)
            optimize!(m)
            set_prop!(mg, vertex, :m, m)
        end
        set_indexing_prop!(mg, :m)
        pdiffs, qdiffs, vdiffs = get_diffs(mg)
        maxp = maximum(abs.(pdiffs))
        maxq = maximum(abs.(qdiffs))
        maxv = maximum(abs.(vdiffs))
        if verbose
            i += 1
            println("\niterate $i")
            println("max pdiff $(maxp)")
            println("max qdiff $(maxq)")
            println("max vdiff $(maxv)")
        end
        # we don't want to set_inputs! if we are about to break
        diff = maximum([maxp, maxq, maxv])
        if diff > abs(tol)
            set_inputs!(mg; α=α)
        end
    end
end


"""
    solve_metagraph!(mg::MetaDiGraph, builder::Function, tols::Vector{T}; α::T=0.5, verbose=false) where T <: Real
    
Given a MetaDiGraph and a JuMP Model `builder` method iteratively solve the models until the `tols` are 
met for the differences provided by `BranchFlowModel.get_diffs`. 

The `builder` must accept only one argument of type `BranchFlowModel.AbstractInputs` that returns 
a `JuMP.AbstractModel`. Each model returned from the `builder` is stored as an `:m` property in 
each vertex of `mg`.

!!! note 
    The `tols` should have a length of three. The first value is compared to the maximum absolute difference
    in Pj, the second for Qj, and the third for |v|. All differences are calculated at the leaf/substation
    connections.
"""
function solve_metagraph!(mg::MetaDiGraph, builder::Function, tols::Vector{R}; α::T=0.5, verbose=false) where {T <: Real, R <: Real}
    init_inputs!(mg)
    tol_not_met = true
    i = 0
    while tol_not_met
        for (p, vertex) in mg[:p]
            m = builder(p)
            optimize!(m)
            set_prop!(mg, vertex, :m, m)
        end
        set_indexing_prop!(mg, :m)
        pdiffs, qdiffs, vdiffs = get_diffs(mg)
        maxp = maximum(abs.(pdiffs))
        maxq = maximum(abs.(qdiffs))
        maxv = maximum(abs.(vdiffs))
        if verbose
            i += 1
            println("\niterate $i")
            println("max pdiff $(maxp)")
            println("max qdiff $(maxq)")
            println("max vdiff $(maxv)")
        end
        if all( val <= tol for (val,tol) in zip([maxp, maxq, maxv], tols))
            tol_not_met = false
            continue
        end
        set_inputs!(mg; α=α)
    end
end


"""
    solve_metagraph!(mg::MetaDiGraph, builder::Dict{Int64, Function}, tols::Vector{T}; α::T=0.5, verbose=false) where T <: Real
    
Given a MetaDiGraph and a JuMP Model `builder` method iteratively solve the models until the `tols` are 
met for the differences provided by `BranchFlowModel.get_diffs`. 
The `builder` dict is used to build each model for the corresponding vertex key.

Each function in the `builder` dict must accept only one argument of type `BranchFlowModel.AbstractInputs` that returns 
a `JuMP.AbstractModel`. Each model returned from the builder function is stored as an `:m` property in 
each vertex of `mg`.

!!! note 
    The `tols` should have a length of three. The first value is compared to the maximum absolute difference
    in Pj, the second for Qj, and the third for |v|. All differences are calculated at the leaf/substation
    connections.
"""
function solve_metagraph!(mg::MetaDiGraph, builder::Dict{Int64, <:Function}, tols::Vector{R}; α::T=0.5, verbose=false) where {T <: Real, R <: Real}
    init_inputs!(mg)
    tol_not_met = true
    i = 0
    while tol_not_met
        for (p, vertex) in mg[:p]
            m = builder[vertex](p)
            optimize!(m)
            set_prop!(mg, vertex, :m, m)
        end
        set_indexing_prop!(mg, :m)
        pdiffs, qdiffs, vdiffs = get_diffs(mg)
        maxp = maximum(abs.(pdiffs))
        maxq = maximum(abs.(qdiffs))
        maxv = maximum(abs.(vdiffs))
        if verbose
            i += 1
            println("\niterate $i")
            println("max pdiff $(maxp)")
            println("max qdiff $(maxq)")
            println("max vdiff $(maxv)")
        end
        if all( val <= tol for (val,tol) in zip([maxp, maxq, maxv], tols))
            tol_not_met = false
            continue
        end
        set_inputs!(mg; α=α)
    end
end
