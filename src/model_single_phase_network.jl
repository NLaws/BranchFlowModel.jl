"""
    build_model!(m::JuMP.AbstractModel, net::Network{SinglePhase})

Add variables and constraints to `m` using the values in `net`. Calls the following functions:
```julia
add_variables(m, net)
constrain_power_balance(m, net)
constrain_substation_voltage(m, net)
constrain_KVL(m, net)
constrain_cone(m, net)
constrain_loads(m, net)
```
"""
function build_model!(m::JuMP.AbstractModel, net::Network{SinglePhase}; relaxed::Bool=true)

    add_variables(m, net)
    constrain_power_balance(m, net)
    constrain_substation_voltage(m, net)
    constrain_KVL(m, net)
    constrain_loads(m, net)
    if relaxed
        constrain_cone(m, net)
    else
        constrain_bilinear(m, net)
    end
end


function add_variables(m, net::Network{SinglePhase})
    bs = collect(busses(net))
    es = collect(edges(net))
    # bus injections are expressions, defined in constrain_power_balance

    # TODO warn when applying power injection lower bounds and SDP, radial b/c for radial
    # networks with power injections unbounded below (voltage angle relaxation) the SOCP (more
    # efficiently) yields the same solution as the SDP and the (non-convex) exact equations.

    # voltage squared
    CommonOPF.add_time_vector_variables!(m, net, :vsqrd, bs)
    for bus in busses(net)
        JuMP.set_lower_bound.(m[:vsqrd][bus], (net.v_lolim)^2)
    end
    # TODO upper bounds
    
    # line flows, net power sent from i to j
    # @variable(m, net.P_lo_bound <= Pij[edges(net), T] <= p.P_up_bound )
    # @variable(m, net.Q_lo_bound <= Qij[edges(net), T] <= p.Q_up_bound )
    CommonOPF.add_time_vector_variables!(m, net, :Pij, es)
    CommonOPF.add_time_vector_variables!(m, net, :Qij, es)

    # current squared (non-negative)
    CommonOPF.add_time_vector_variables!(m, net, :lij, es)

    # @constraint(m, [edge in net.edge_keys, t in T],
    #     lij[edge, t] <= net.amps_limit[edge]
    # )

    nothing
end


"""
    function constrain_power_balance(m, net::Network)

Pij in - losses == sum of line flows out + net injection
NOTE: using sum over Pij for future expansion to mesh grids
i -> j -> k
"""
function constrain_power_balance(m, net::Network)
    Pij = m[:Pij]
    Qij = m[:Qij]
    lij = m[:lij]
    Pj = m[:Pj] = Dict()  # we fill these in with expression indexed on bus and time
    Qj = m[:Qj] = Dict()
    # NOTE with Pj and Qj as expressions (instead of variables) the # of variables is reduced by
    # (Nnodes - 1)*8760 and number of constraints by 6*(Nnodes - 1)*8760 
    
    # TODO using expressions for net injections means that we do not have to delete and redifine
    # constraints? i.e. can just define another expression for Pj when it is a function of a DER ?
    for j in busses(net)

        # source nodes, injection = flows out
        if isempty(i_to_j(j, net)) && !isempty(j_to_k(j, net))
            substation_Pload = zeros(net.Ntimesteps)
            if j in real_load_busses(net)
                substation_Pload = net[j][:Load][:kws1] * 1e3 / net.Sbase
            end
            Pj[j] = @expression(m, [t = 1:net.Ntimesteps],
                - substation_Pload[t] - sum( Pij[(j,k)][t] for k in j_to_k(j, net) )
            )
            substation_Qload = zeros(net.Ntimesteps)
            if j in reactive_load_busses(net)
                substation_Qload = net[j][:Load][:kvars1] * 1e3 / net.Sbase
            end
            Qj[j] = @expression(m, [t = 1:net.Ntimesteps],
                - substation_Qload[t] - sum( Qij[(j,k)][t] for k in j_to_k(j, net) )
            )
        
        # unconnected nodes
        elseif isempty(i_to_j(j, net)) && isempty(j_to_k(j, net))
            @warn "Bus $j has no edges, setting Pj and Qj to zero."
            Pj[j] = zeros(net.Ntimesteps)
            Qj[j] = zeros(net.Ntimesteps)

        # leaf nodes / sinks, flows in = draw out
        elseif !isempty(i_to_j(j, net)) && isempty(j_to_k(j, net))
            Pj[j] = @expression(m, [t = 1:net.Ntimesteps],
                -sum( Pij[(i,j)][t] for i in i_to_j(j, net) ) +
                sum( lij[(i,j)][t] * rij(i,j,net) for i in i_to_j(j, net) ) 
            )
            Qj[j] = @expression(m, [t = 1:net.Ntimesteps],
                -sum( Qij[(i,j)][t] for i in i_to_j(j, net) ) +
                sum( lij[(i,j)][t] * xij(i,j,net) for i in i_to_j(j, net) )
                # + p.shunt_susceptance[j] * m[:vsqrd][j,t]
            )
        
        # intermediate nodes
        else
            Pj[j] = @expression(m, [t = 1:net.Ntimesteps],
                -sum( Pij[(i,j)][t] for i in i_to_j(j, net) )
                + sum( lij[(i,j)][t] * rij(i,j,net) for i in i_to_j(j, net) ) 
                + sum( Pij[(j,k)][t] for k in j_to_k(j, net) )
            )
            Qj[j] = @expression(m, [t = 1:net.Ntimesteps],
                -sum( Qij[(i,j)][t] for i in i_to_j(j, net) ) 
                + sum( lij[(i,j)][t] * xij(i,j,net) for i in i_to_j(j, net) )
                + sum( Qij[(j,k)][t] for k in j_to_k(j, net) ) 
                # + p.shunt_susceptance[j] * m[:vsqrd][j,t]
            )
        end
    end
    # TODO Farivar and Low have b*v and q*v in these equations for shunts? neglected for now

    nothing
end


function constrain_substation_voltage(m, net::Network)
    if typeof(net.v0) <: Real
        @constraint(m, con_substationV[t = 1:net.Ntimesteps],
            m[:vsqrd][net.substation_bus][t] == net.v0^2
        )
    else  # vector of time
        @constraint(m, con_substationV[t = 1:net.Ntimesteps],
            m[:vsqrd][net.substation_bus][t] == net.v0[t]^2
        )
    end
    nothing
end


function constrain_KVL(m, net::Network)
    w = m[:vsqrd]
    P = m[:Pij]
    Q = m[:Qij]
    l = m[:lij]
    m[:vcons] = Dict()
    for j in busses(net)
        for i in i_to_j(j, net)  # for radial network there is only one i in i_to_j
            # if !( (i,j) in keys(p.regulators) )
                rᵢⱼ = rij(i,j,net)
                xᵢⱼ = xij(i,j,net)
                m[:vcons][j] = @constraint(m, [t = 1:net.Ntimesteps],
                    w[j][t] == w[i][t]
                        - 2*(rᵢⱼ * P[(i,j)][t] + xᵢⱼ * Q[(i,j)][t])
                        + (rᵢⱼ^2 + xᵢⱼ^2) * l[(i,j)][t]
                )
            # else
            #     if has_vreg(p, j)
            #         m[:vcons][j] = @constraint(m, [t = 1:net.Ntimesteps],
            #             w[j,t] == p.regulators[(i,j)][:vreg]^2
            #         )
            #     else  # default turn_ratio is 1.0
            #         m[:vcons][j] = @constraint(m, [t = 1:net.Ntimesteps],
            #             w[j,t] == w[i,t] * p.regulators[(i,j)][:turn_ratio]^2 
            #         )
            #     end
            # end
        end
    end
    nothing
end


function constrain_cone(m, net::Network)
    w = m[:vsqrd]
    P = m[:Pij]
    Q = m[:Qij]
    l = m[:lij]
    for j in busses(net)
        for i in i_to_j(j, net)  # for radial network there is only one i in i_to_j
            # # equivalent but maybe not handled as well in JuMP ?
            # @constraint(m, [t = 1:net.Ntimesteps],
            #     w[i,t] * l[(i,j), t] ≥ P[(i,j),t]^2 + Q[(i,j),t]^2
            # )
            @constraint(m, [t = 1:net.Ntimesteps], 
                [w[i][t]/2, l[(i,j)][t], P[(i,j)][t], Q[(i,j)][t]] in JuMP.RotatedSecondOrderCone()
            )
        end
    end
end


function constrain_bilinear(m, net::Network)
    w = m[:vsqrd]
    P = m[:Pij]
    Q = m[:Qij]
    l = m[:lij]
    for j in busses(net)
        for i in i_to_j(j, net)  # for radial network there is only one i in i_to_j
            @constraint(m, [t = 1:net.Ntimesteps],
                w[i][t] * l[(i,j)][t] == P[(i,j)][t]^2 + Q[(i,j)][t]^2
            )
        end
    end
end


"""
    constrain_loads(m, net::Network)

- set net injections Pj/Qj to negative of Inputs.Pload/Qload, which are normalized by Sbase when creating Inputs
- keys of P/Qload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus

TODO this method is same as LinDistFlow single phase: should it be moved to CommonOPF?
"""
function constrain_loads(m, net::Network)
    Pj = m[:Pj]
    Qj = m[:Qj]
    m[:injectioncons] = Dict()
    for j in setdiff(busses(net), [net.substation_bus])
        m[:injectioncons][j] = Dict()
        if j in real_load_busses(net)
            con = @constraint(m, [t = 1:net.Ntimesteps],
                Pj[j][t] == -net[j][:Load][:kws1][t] * 1e3 / net.Sbase
            )
        else
            con = @constraint(m, [t = 1:net.Ntimesteps],
                Pj[j][t] == 0
            )
        end
        m[:injectioncons][j]["p"] = con
        if j in reactive_load_busses(net)
            con = @constraint(m, [t = 1:net.Ntimesteps],
                Qj[j][t] == -net[j][:Load][:kvars1][t] * 1e3 / net.Sbase
            )
        else
            con = @constraint(m, [t = 1:net.Ntimesteps],
                Qj[j][t] == 0
            )
        end
        m[:injectioncons][j]["q"] = con
    end
    nothing
end
