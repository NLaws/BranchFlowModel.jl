



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


function get_bus_values(var::Symbol, m::JuMP.AbstractModel, p::Inputs{SinglePhase})
    vals = value.(m[var])
    d = Dict()
    for b in p.busses
        d[b] = vals[b,:].data
        if var == :vsqrd || var == :lij
            d[b] = sqrt.(d[b])
        end
    end
    return d
end


function metagraph_voltages(mg::MetaDiGraph)
    d = Dict()
    for v in vertices(mg)
        merge!(d, get_bus_values(:vsqrd, mg[v, :m], mg[v, :p]))
    end
    return d
end