
function get_constraints_by_variable_name(m::JuMP.AbstractModel, v::AbstractString)
    ac = ConstraintRef[]
    for tup in list_of_constraint_types(m)
        append!(ac, all_constraints(m, tup[1], tup[2]))
    end
    filter( cr -> occursin(v, string(cr)), ac )
end


"""
    injection_at_bus(j::String, net::Network{MultiPhase})::Tuple{Vector{Vector{<:Real}}, Vector{Vector{<:Real}}}

return the real and reactive power injections as vectors with 3 phase indices and net.Ntimesteps time
indices like:
```julia
Pj, Qj = injection_at_bus(my_bus, net)
...
Pj[phase][time_step]
```
"""
function injection_at_bus(j::String, net::Network{MultiPhase})::Tuple{Vector{Vector{<:Real}}, Vector{Vector{<:Real}}}
    Pj = [zeros(net.Ntimesteps) for _ in 1:3] # first dim is phase, like Pj[phs][t]
    Qj = [zeros(net.Ntimesteps) for _ in 1:3]
    if j in real_load_busses(net)
        for phs in 1:3
            Pj[phs] = -net[j, :kws, phs] * 1e3 / net.Sbase
        end
    end
    if j in reactive_load_busses(net)
        for phs in 1:3 
            Qj[phs] = -net[j, :kvars, phs] * 1e3 / net.Sbase
        end
    end
    return Pj, Qj
end
