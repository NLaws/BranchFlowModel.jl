abstract type AbstractResults end

"""
TODO
"""
struct Results <: AbstractResults
    voltage_magnitudes
    real_power_injections
    reactive_power_injections
    current_magnitudes
    real_sending_end_powers
    reactive_sending_end_powers
    shadow_prices
end


"""
    Results(m::AbstractModel, p::Inputs{SinglePhase}; digits=8)

return a `Results` struct with fieldnames:

    voltage_magnitudes
    real_power_injections
    reactive_power_injections
    current_magnitudes
    real_sending_end_powers
    reactive_sending_end_powers

"""
function Results(m::AbstractModel, p::Inputs{SinglePhase}; digits=8)

    vs  = get_variable_values(:vsqrd, m, p, digits=digits)
    Pj  = get_variable_values(:Pj,    m, p, digits=digits)
    Qj  = get_variable_values(:Qj,    m, p, digits=digits)
    lij = get_variable_values(:lij,   m, p, digits=digits)
    Pij = get_variable_values(:Pij,   m, p, digits=digits)
    Qij = get_variable_values(:Qij,   m, p, digits=digits)
    prices = Dict(
        j => JuMP.dual.(m[:loadbalcons][j]["p"])
        for j in p.busses
    )

    Results(vs, Pj, Qj, lij, Pij, Qij, prices)
end


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


function metagraph_voltages(mg::MetaDiGraph)
    d = Dict()
    for v in vertices(mg)
        merge!(d, get_variable_values(:vsqrd, mg[v, :m], mg[v, :p]))
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
