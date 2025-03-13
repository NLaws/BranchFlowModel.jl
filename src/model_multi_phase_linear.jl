
# TODO mv all linear stuff to this file and include the file
# TODO MP/Qij should have linear in their names
"""
    MPij(i::String, j::String, net::Network{MultiPhase})

Real power coefficients for 3 phase voltage drop from node i to j
"""
function MPij(i::String, j::String, net::Network{MultiPhase})
    M = zeros((3,3))
    r = rij(i,j,net)
    x = xij(i,j,net)
    M[1,:] = [-2r[1,1]             r[1,2]-sqrt(3)x[1,2] r[1,3]+sqrt(3)x[1,3]]
    M[2,:] = [r[2,1]+sqrt(3)x[2,1] -2r[2,2]             r[2,3]-sqrt(3)x[2,3]]
    M[3,:] = [r[3,1]-sqrt(3)x[3,1] r[3,2]+sqrt(3)x[3,2] -2r[3,3]            ]
    return M
end


"""
    MQij(i::String, j::String, net::Network{MultiPhase})

Reactive power coefficients for 3 phase voltage drop from node i to j
"""
function MQij(i::String, j::String, net::Network{MultiPhase})
    M = zeros((3,3))
    r = rij(i,j,net)
    x = xij(i,j,net)
    M[1,:] = [-2x[1,1]             x[1,2]+sqrt(3)r[1,2] x[1,3]-sqrt(3)r[1,3]]
    M[2,:] = [x[2,1]-sqrt(3)r[2,1] -2x[2,2]             x[2,3]+sqrt(3)r[2,3]]
    M[3,:] = [x[3,1]+sqrt(3)r[3,1] x[3,2]-sqrt(3)r[3,2] -2x[3,3]            ]
    return M
end


function constrain_KVL_linear(m, net::Network{MultiPhase})
    w = m[:vsqrd]
    p = m[:pij]
    q = m[:qij]
    for j in busses(net)  # substation_bus in here but has empty i_to_j(j, net)
        for i in i_to_j(j, net)  # for radial network there is only one i in i_to_j
            i_j = (i, j)

            if !( isa(net[(i,j)], CommonOPF.VoltageRegulator) )
                
                MP = MPij(i,j,net)
                MQ = MQij(i,j,net)

                for phs in phases_into_bus(net, j)
                    @constraint(m, [t in 1:net.Ntimesteps],
                        w[j][t][phs] == w[i][t][phs]
                            + sum(MP[phs,k] * p[i_j][t][k] for k=phases_into_bus(net, j))
                            + sum(MQ[phs,k] * q[i_j][t][k] for k=phases_into_bus(net, j))
                    )
                end

            else
                reg = net[(i,j)]
                if !ismissing(reg.vreg_pu)
                    @constraint(m, [t in 1:net.Ntimesteps],
                        w[j][t] .== reg.vreg_pu.^2
                    )
                else  # use turn ratio
                    @constraint(m, [t in 1:net.Ntimesteps],
                        w[j][t] .== w[i][t] * reg.turn_ratio^2 
                    )
                end
            end
        end
    end
    nothing
end
