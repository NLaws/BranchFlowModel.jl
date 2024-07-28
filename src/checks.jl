# """
#     check_soc_inequalities(m::JuMP.AbstractModel, p::Inputs)

# create and return a vector of the gaps in the second order cone constraints
# """
# function check_soc_inequalities(m::JuMP.AbstractModel, p::Inputs)
#     @assert termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]
#     w = value.(m[:vsqrd])
#     lij = value.(m[:lij])
#     Pj = value.(m[:Pj])
#     Qj = value.(m[:Qj])
#     gaps = Vector{Float64}()  # todo size empty array based on number of SOC cons

#     # TODO mesh accounting
#     for j in p.busses, i in i_to_j(j, p), t in 1:p.Ntimesteps
#         push!(gaps, 
#             lij[string(i*"-"*j), t] * w[i, t] - Pj[j,t]^2 + Qj[j,t]^2 
#         )
#     end
#     return gaps
# end


"""
    check_rank_one(m::JuMP.AbstractModel, net::Network, tol=1e-3)

Check the rank of the `m[:H]` matrices from the PSD cone constraints.
Warnings express any values with rank greater than one.
"""
function check_rank_one(m::JuMP.AbstractModel, net::Network, tol=1e-3)
    for j in busses(net), t in 1:net.Ntimesteps
        if j == net.substation_bus continue end
        # TODO use LinearAlgebra.rank w/ atol kwarg
        eigs = eigvals!(JuMP.value.(m[:H][t][j]))
        eigs ./= maximum(eigs)
        rank_H = sum(map(x-> x > tol ? 1 : 0, eigs))
        if rank_H > 1
            @warn("Bus $j in time step $t has H matrix of rank $rank_H")
        end
    end
end
