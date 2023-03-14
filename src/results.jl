



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
