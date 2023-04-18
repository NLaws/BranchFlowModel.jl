



function voltage_values_by_time_bus(m::JuMP.AbstractModel, p::Inputs{BranchFlowModel.MultiPhase})
    d = Dict{Int, Dict{String, AbstractMatrix}}()
    for t in 1:p.Ntimesteps
        d[t] = Dict()
        for b in p.busses
            d[t][b] = sqrt.(value.(m[:w][t][b]))
        end
    end
    return d
end



function current_values_by_time_edge(m::JuMP.AbstractModel, p::Inputs{BranchFlowModel.MultiPhase})
    d = Dict{Int, Dict{String, AbstractMatrix}}()
    for t in 1:p.Ntimesteps
        d[t] = Dict()
        for e in p.edge_keys
            d[t][e] = sqrt.(value.(m[:l][t][e]))
        end
    end
    return d
end


function line_flow_values_by_time_edge(m::JuMP.AbstractModel, p::Inputs{BranchFlowModel.MultiPhase})
    d = Dict{Int, Dict{String, AbstractMatrix}}()
    for t in 1:p.Ntimesteps
        d[t] = Dict()
        for e in p.edge_keys
            d[t][e] = value.(m[:Sij][t][e])
        end
    end
    return d
end


function get_variable_values(var::Symbol, m::JuMP.AbstractModel, p::Inputs{SinglePhase})
    d = Dict()
    if var in [:Pj, :Qj, :vsqrd]  # TODO make these a const in CommonOPF
        vals = value.(m[var])
        for b in p.busses
            d[b] = vals[b,:].data
            if var == :vsqrd
                d[b] = sqrt.(d[b])
            end
        end
    elseif var in [:Pij, :Qij, :lij]  # TODO make these a const in CommonOPF
        vals = value.(m[var])
        for ek in p.edge_keys
            d[ek] = vals[ek,:].data
            if var == :lij
                d[ek] = sqrt.(d[ek])
            end
        end
    else
        @warn "$var is not a valid variable symbol"
    end
    return d
end


function get_bus_values(var::Symbol, m::JuMP.AbstractModel, p::Inputs{SinglePhase})
    @warn "get_bus_values will be deprecated in favor of get_variable_values in the next major release."
    get_variable_values(var, m, p)
end


function metagraph_voltages(mg::MetaDiGraph)
    d = Dict()
    for v in vertices(mg)
        merge!(d, get_bus_values(:vsqrd, mg[v, :m], mg[v, :p]))
    end
    return d
end


function get_edge_values(var_prefix::AbstractString, m::JuMP.AbstractModel, p::Inputs)
    vals = Float64[]
    for edge in p.edges
        var = string(var_prefix, "[", edge[1], "-", edge[2], "]")
        try
            val = value(variable_by_name(m, var))
            if startswith(var_prefix, "l")
                val = sqrt(val)
            end
            push!(vals, round(val; digits=8))
        catch e
            println(var, "failed", e)
        end
    end
    return vals
end
