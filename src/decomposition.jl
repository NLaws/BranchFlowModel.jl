
"""
    connecting_busses(mg::MetaGraphsNext.MetaGraph, v)

return Dict with keys for out-vertices of node `v` and values for their substation_bus
"""
function connecting_busses(mg::MetaGraphsNext.MetaGraph, v)
    cbusses = Dict()
    for neighb in  outneighbors(mg, v)
        if mg[neighb, :p].substation_bus in mg[v, :p].busses
            cbusses[neighb] = mg[neighb, :p].substation_bus
        end
    end
    return cbusses
end


"""
    leaf_vertices(mg::MetaGraphsNext.MetaGraph)

returns `Vector{Int64}` containing all of the leaf vertices in `mg`
"""
function leaf_vertices(mg::MetaGraphsNext.MetaGraph)
    leafs = Int64[]
    for v in vertices(mg)
        if outdegree(mg, v) == 0
            push!(leafs, v)
        end
    end
    return leafs
end


"""
    set_inputs!(mg::MetaGraphsNext.MetaGraph; α::Float64=0.0)

Set the shared values in each subgraph / vertex of mg:
1. set the current vertex's v0 to its inneighbor's voltage
2. set the current vertex P/Qload to the outneighbors' substation_bus loads
"""
function set_inputs!(mg::MetaGraphsNext.MetaGraph; α::R=0.0) where R <: Real
    for v in mg.graph_data[:load_sum_order] # ~breadth first search of vertices
        # if v has inneighbors use their voltages at connections
        p_below = mg[v]
        for v_above in inneighbors(mg, v)
            m_above = get_prop(mg, v_above, :m)
            if p_below.substation_bus in reg_busses(p_below) || α == 0.0 
                # then p_below.v0 is exactly the voltage at the same bus in graph above 
                p_below.v0 = sqrt.(value.(m_above[:vsqrd][p_below.substation_bus, :])).data
            else # use weighted average of new and old values
                m_below = get_prop(mg, v, :m)
                v_kp1 = sqrt.(value.(m_above[:vsqrd][p_below.substation_bus, :])).data
                v_k   = sqrt.(value.(m_below[:vsqrd][p_below.substation_bus, :])).data
                p_below.v0 = (v_kp1 + α .* v_k) ./ (1 + α)
            end
        end
        # if v has outneighbors then use v_below's m[:Pj][p_below.substation_bus,:] * p_above.Sbase as v's Pload at the same bus
        # have set the loads from deepest vertices to shallowest to get correct sums
        p_above = mg[v]
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
    breadth_first_nodes_from_leafs(mg:::MetaGraphsNext.MetaGraph)

returns a `Vector{Int64}` for the vertices in `mg` from the leafs to the source.
"""
function breadth_first_nodes_from_leafs(mg::MetaGraphsNext.MetaGraph)
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


function check_statuses(mg::MetaGraphsNext.MetaGraph)
    good_status = [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]
    for v in vertices(mg)
        if !(termination_status(mg[v, :m]) in good_status)
            @warn("vertex $v status: $(termination_status(mg[v, :m]))")
        end
    end
end



"""
    get_diffs(mg::MetaGraphsNext.MetaGraph)

Uses the JuMP Models stored in mg[:m] to calculate the difference between Pj, Qj, and |v| at every
leaf/substation connection. 

returns three `Float64[]`
"""
function get_diffs(mg::MetaGraphsNext.MetaGraph)
    pdiffs, qdiffs, vdiffs = Float64[], Float64[], Float64[]
    for v in mg.graph_data[:load_sum_order] # ~breadth first search of vertices
        # if v has inneighbors use their voltages at connections
        p_below = mg[v]
        for v_above in inneighbors(mg, v)
            m_above = get_prop(mg, v_above, :m)
            volts_above = sqrt.(value.(m_above[:vsqrd][p_below.substation_bus, :]))  # vector of time
            if has_vreg(p_below, p_below.substation_bus)
                # the voltage at p_below.substation_bus is fixed and known, essentially starts new problem
                # however, the NLP can lead to violation of the substation voltage constraint
                v_below = sqrt.(value.(mg[v, :m][:vsqrd][p_below.substation_bus, :])).data
                push!(vdiffs, sum(vreg(p_below, p_below.substation_bus) .- v_below) / p_below.Ntimesteps)
            else
                push!(vdiffs, sum(abs.(p_below.v0 .- volts_above)) / p_below.Ntimesteps)
            end
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
    metagraph_to_json(mg::MetaGraphsNext.MetaGraph, filename::String)

Dump the vertices, edges, and the busses and lines for each vertex in the metagraph to JSON
"""
function metagraph_to_json(mg::MetaGraphsNext.MetaGraph, filename::String)
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
    solve_metagraph!(mg::MetaGraphsNext.MetaGraph, builder::Function, tol::T; α::T=0.5, verbose=false) where T <: Real

Given a MetaGraphsNext.MetaGraph and a JuMP Model `builder` method iteratively solve the models until the `tol` is 
met for the differences provided by `BranchFlowModel.get_diffs`. 

The `builder` must accept only one argument of type `BranchFlowModel.AbstractInputs` that returns 
a `JuMP.AbstractModel`. Each model returned from the `builder` is stored as an `:m` property in 
each vertex of `mg`.

!!! note 
    `tol` is compared to the maximum absolute value of all the p, q, and v differences.
"""
function solve_metagraph!(mg::MetaGraphsNext.MetaGraph, builder::Function, tol::R; α::T=0.5, verbose=false) where {T <: Real, R <: Real}
    CommonOPF.init_split_networks!(mg)
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
            println()
            check_statuses(mg)
            println("iterate $i")
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
    solve_metagraph!(mg::MetaGraphsNext.MetaGraph, builder::Function, tols::Vector{T}; α::T=0.5, verbose=false) where T <: Real
    
Given a MetaGraphsNext.MetaGraph and a JuMP Model `builder` method iteratively solve the models until the `tols` are 
met for the differences provided by `BranchFlowModel.get_diffs`. 

The `builder` must accept only one argument of type `BranchFlowModel.AbstractInputs` that returns 
a `JuMP.AbstractModel`. Each model returned from the `builder` is stored as an `:m` property in 
each vertex of `mg`.

!!! note 
    The `tols` should have a length of three. The first value is compared to the maximum absolute difference
    in Pj, the second for Qj, and the third for |v|. All differences are calculated at the leaf/substation
    connections.
"""
function solve_metagraph!(mg::MetaGraphsNext.MetaGraph, builder::Function, tols::Vector{R}; α::T=0.5, verbose=false) where {T <: Real, R <: Real}
    CommonOPF.init_split_networks!(mg)
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
            println()
            check_statuses(mg)
            println("iterate $i")
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
    solve_metagraph!(mg::MetaGraphsNext.MetaGraph, builder::Dict{Int64, Function}, tols::Vector{T}; α::T=0.5, verbose=false) where T <: Real
    
Given a MetaGraphsNext.MetaGraph and a JuMP Model `builder` method iteratively solve the models until the `tols` are 
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
function solve_metagraph!(mg::MetaGraphsNext.MetaGraph, builder::Dict{Int64, <:Function}, tols::Vector{R}; α::T=0.5, verbose=false) where {T <: Real, R <: Real}
    CommonOPF.init_split_networks!(mg)
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
            println()
            check_statuses(mg)
            println("iterate $i")
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
