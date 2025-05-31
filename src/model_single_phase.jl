"""
    build_bfm!(m::JuMP.AbstractModel, net::Network{SinglePhase}, mtype::ModelType=Unrelaxed)

Top-level single phase builder that dispatches the ModelType enum
"""
function build_bfm!(m::JuMP.AbstractModel, net::Network{SinglePhase}, mtype::ModelType=Unrelaxed)
    build_bfm!(m::JuMP.AbstractModel, net::Network{SinglePhase}, Val(mtype))
end


"""
    build_bfm!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{Unrelaxed})

Add variables and constraints to `m` using the values in `net`. Calls the following functions:
```julia
add_linear_variables(m, net)
add_vsqrd_variables(m, net)
add_isqrd_variables(m, net)
constrain_power_balance_with_isqrd_losses(m, net)
constrain_substation_voltage(m, net)
constrain_KVL(m, net)
constrain_bilinear(m, net)
```
"""
function build_bfm!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{Unrelaxed})
    # TODO the Unrelaxed model has the angle relaxation, should have both unrelaxed and with angle relaxation?
    add_linear_variables(m, net)
    add_vsqrd_variables(m, net)
    add_isqrd_variables(m, net)
    constrain_power_balance_with_isqrd_losses(m, net)
    constrain_substation_voltage(m, net)
    constrain_KVL(m, net)
    constrain_bilinear(m, net)
end


"""
    build_bfm!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{SecondOrderCone})

Add variables and constraints to `m` using the values in `net`. Calls the following functions:
```julia
add_linear_variables(m, net)
add_vsqrd_variables(m, net)
add_isqrd_variables(m, net)
constrain_power_balance_with_isqrd_losses(m, net)
constrain_substation_voltage(m, net)
constrain_KVL(m, net)
constrain_cone(m, net)
```
"""
function build_bfm!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{SecondOrderCone})
    add_linear_variables(m, net)
    add_vsqrd_variables(m, net)
    add_isqrd_variables(m, net)
    constrain_power_balance_with_isqrd_losses(m, net)
    constrain_substation_voltage(m, net)
    constrain_KVL(m, net)
    constrain_cone(m, net)
end


"""
    build_bfm!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{Linear})

Add variables and constraints to `m` using the values in `net`. Calls the following functions:
```julia
add_linear_variables(m, net)
add_vsqrd_variables(m, net)
add_isqrd_variables(m, net)
constrain_power_balance_with_isqrd_losses(m, net)
constrain_substation_voltage(m, net)
constrain_KVL(m, net)
constrain_cone(m, net)
```
"""
function build_bfm!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{Linear})
    add_linear_variables(m, net)
    add_vsqrd_variables(m, net)
    constrain_power_balance_linear(m, net)
    constrain_substation_voltage(m, net)
    constrain_KVL_linear(m, net)
end


"""
    add_linear_variables(m, net::Network{SinglePhase})

Add variables for the single-phase, linear model:
- `pij` and `qij` for all `edges(net)`
- `p0` and `q0` slack bus power
"""
function add_linear_variables(m, net::Network{SinglePhase})
    es = collect(edges(net))
    
    # line flows, net power sent from i to j
    CommonOPF.add_time_vector_variables!(m, net, :pij, es)
    CommonOPF.add_time_vector_variables!(m, net, :qij, es)

    # TODO slack bus variables in CommonOPF
    # slack bus variables
    @variable(m, p0[1:net.Ntimesteps])
    @variable(m, q0[1:net.Ntimesteps])

    nothing
end


function add_vsqrd_variables(m, net::Network{SinglePhase})
    bs = collect(busses(net))

    # TODO warn when applying power injection lower bounds and SDP, radial b/c for radial
    # networks with power injections unbounded below (voltage angle relaxation) the SOCP (more
    # efficiently) yields the same solution as the SDP and the (non-convex) exact equations.

    # voltage squared
    CommonOPF.add_time_vector_variables!(m, net, :vsqrd, bs)
    for bus in busses(net)
        if !ismissing(net.bounds.v_lower_mag)
            JuMP.set_lower_bound.(m[:vsqrd][bus], (net.bounds.v_lower_mag)^2)
        end
        if !ismissing(net.bounds.v_upper_mag)
            JuMP.set_upper_bound.(m[:vsqrd][bus], (net.bounds.v_upper_mag)^2)
        end
    end
    # TODO more bounds in a centralized fashion

    nothing
end


function add_isqrd_variables(m, net::Network{SinglePhase})
    es = collect(edges(net))

    # current squared (non-negative)
    CommonOPF.add_time_vector_variables!(m, net, :lij, es)
    for edge in es
        JuMP.set_lower_bound.(m[:lij][edge], 0.0)
    end

    # @constraint(m, [edge in net.edge_keys, t in T],
    #     lij[edge, t] <= net.amps_limit[edge]
    # )

    nothing

end


"""
    function constrain_power_balance_with_isqrd_losses(m, net::Network)

Define the m[:loadbalcons][bus] ∀ bus ∈ busses(net) as a Dict of constraints. The keys are "p" and
"q" for real and reactive power balance respectively. The values are the JuMP constraints.

    ∑ pij in - losses + net injection - ∑ Pjk out = 0

The net injection are user defined loads. If one wishes to make the net injection a decision
variable then delete the constraint and redefine the constraint with your decision variable.

NOTE: using sum over pij for future expansion to mesh grids and the convention:
i -> j -> k
"""
function constrain_power_balance_with_isqrd_losses(m, net::Network{SinglePhase})
    pij = m[:pij]
    qij = m[:qij]
    lij = m[:lij]

    m[:loadbalcons] = Dict()

    m[:loadbalcons][net.substation_bus] = Dict()

    # NOTE defining flows in as positive and flows out as negative makes the marginal prices positive
    
    for j in busses(net)

        m[:loadbalcons][j] = Dict()

        # check for shunt admittance
        shunt_susceptance = 0.0
        if :ShuntAdmittance in keys(net[j])
            # flip sign of susceptance b/c we want complex conjugate
            shunt_susceptance = -net[j][:ShuntAdmittance].b
        end

        # check for capacitors (only fixed for now)
        var = 0.0
        if :Capacitor in keys(net[j])
            var = net[j][:Capacitor].var
        end

        # check for loads
        pj = zeros(net.Ntimesteps)
        qj = zeros(net.Ntimesteps)
        if j in real_load_busses(net)
            pj = -net[j][:Load].kws1 * 1e3 / net.Sbase
        end
        if j in reactive_load_busses(net)
            qj = -net[j][:Load].kvars1 * 1e3 / net.Sbase
        end

        # define the constraints
        # source nodes, injection = flows out
        if isempty(i_to_j(j, net)) && !isempty(j_to_k(j, net))

            if j == net.substation_bus  # include the slack power variables
                m[:loadbalcons][j]["p"] = @constraint(m, [t = 1:net.Ntimesteps],
                    m[:p0][t] + pj[t] - sum( pij[(j,k)][t] for k in j_to_k(j, net) ) == 0
                )
                m[:loadbalcons][j]["q"] = @constraint(m, [t = 1:net.Ntimesteps],
                    m[:q0][t] + qj[t] - sum( qij[(j,k)][t] for k in j_to_k(j, net) ) == 0 
                )
            else  # a source node with known injection
                m[:loadbalcons][j]["p"] = @constraint(m, [t = 1:net.Ntimesteps],
                    pj[t] - sum( pij[(j,k)][t] for k in j_to_k(j, net) ) == 0
                )
                m[:loadbalcons][j]["q"] = @constraint(m, [t = 1:net.Ntimesteps],
                    qj[t] - sum( qij[(j,k)][t] for k in j_to_k(j, net) ) == 0
                )
            end
        
        # unconnected nodes (should not exist so we warn)
        elseif isempty(i_to_j(j, net)) && isempty(j_to_k(j, net))
            @warn "Bus $j has no edges, setting loadbalcons to zero."
            m[:loadbalcons][j]["p"] = zeros(net.Ntimesteps)
            m[:loadbalcons][j]["q"] = zeros(net.Ntimesteps)

        # leaf nodes / sinks, flows in = draw out
        elseif !isempty(i_to_j(j, net)) && isempty(j_to_k(j, net))
            m[:loadbalcons][j]["p"] = @constraint(m, [t = 1:net.Ntimesteps],
                sum( pij[(i,j)][t] for i in i_to_j(j, net) )
                - sum( lij[(i,j)][t] * rij_per_unit(i,j,net) for i in i_to_j(j, net) ) 
                + pj[t] == 0
            )
            m[:loadbalcons][j]["q"] = @constraint(m, [t = 1:net.Ntimesteps],
                sum( qij[(i,j)][t] for i in i_to_j(j, net) )
                - sum( lij[(i,j)][t] * xij_per_unit(i,j,net) for i in i_to_j(j, net) )
                - shunt_susceptance * m[:vsqrd][j][t]
                + qj[t] + var == 0
            )
        
        # intermediate nodes
        else
            m[:loadbalcons][j]["p"] = @constraint(m, [t = 1:net.Ntimesteps],
                sum( pij[(i,j)][t] for i in i_to_j(j, net) )
                - sum( lij[(i,j)][t] * rij_per_unit(i,j,net) for i in i_to_j(j, net) ) 
                - sum( pij[(j,k)][t] for k in j_to_k(j, net) )
                + pj[t] == 0
            )
            m[:loadbalcons][j]["q"] = @constraint(m, [t = 1:net.Ntimesteps],
                sum( qij[(i,j)][t] for i in i_to_j(j, net) ) 
                - sum( lij[(i,j)][t] * xij_per_unit(i,j,net) for i in i_to_j(j, net) )
                - sum( qij[(j,k)][t] for k in j_to_k(j, net) ) 
                - shunt_susceptance * m[:vsqrd][j][t]
                + qj[t] + var == 0
            )
        end
    end

    nothing
end


"""
    function constrain_power_balance_linear(m, net::Network)

Define the m[:loadbalcons][bus] ∀ bus ∈ busses(net) as a Dict of constraints. The keys are "p" and
"q" for real and reactive power balance respectively. The values are the JuMP constraints.

    ∑ pij in - losses + net injection - ∑ Pjk out = 0

The net injection are user defined loads. If one wishes to make the net injection a decision
variable then delete the constraint and redefine the constraint with your decision variable.

NOTE: using sum over pij for future expansion to mesh grids and the convention:
i -> j -> k
"""
function constrain_power_balance_linear(m, net::Network{SinglePhase})
    pij = m[:pij]
    qij = m[:qij]

    m[:loadbalcons] = Dict()

    m[:loadbalcons][net.substation_bus] = Dict()

    # NOTE defining flows in as positive and flows out as negative makes the marginal prices positive
    
    for j in busses(net)

        m[:loadbalcons][j] = Dict()

        # check for shunt admittance
        shunt_susceptance = 0.0
        if :ShuntAdmittance in keys(net[j])
            # flip sign of susceptance b/c we want complex conjugate
            shunt_susceptance = -net[j][:ShuntAdmittance].b
        end

        # check for capacitors (only fixed for now)
        var = 0.0
        if :Capacitor in keys(net[j])
            var = net[j][:Capacitor].var
        end

        # check for loads
        pj = zeros(net.Ntimesteps)
        qj = zeros(net.Ntimesteps)
        if j in real_load_busses(net)
            pj = -net[j][:Load].kws1 * 1e3 / net.Sbase
        end
        if j in reactive_load_busses(net)
            qj = -net[j][:Load].kvars1 * 1e3 / net.Sbase
        end

        # define the constraints
        # source nodes, injection = flows out
        if isempty(i_to_j(j, net)) && !isempty(j_to_k(j, net))

            if j == net.substation_bus  # include the slack power variables
                m[:loadbalcons][j]["p"] = @constraint(m, [t = 1:net.Ntimesteps],
                    m[:p0][t] + pj[t] - sum( pij[(j,k)][t] for k in j_to_k(j, net) ) == 0
                )
                m[:loadbalcons][j]["q"] = @constraint(m, [t = 1:net.Ntimesteps],
                    m[:q0][t] + qj[t] - sum( qij[(j,k)][t] for k in j_to_k(j, net) ) == 0 
                )
            else  # a source node with known injection
                m[:loadbalcons][j]["p"] = @constraint(m, [t = 1:net.Ntimesteps],
                    pj[t] - sum( pij[(j,k)][t] for k in j_to_k(j, net) ) == 0
                )
                m[:loadbalcons][j]["q"] = @constraint(m, [t = 1:net.Ntimesteps],
                    qj[t] - sum( qij[(j,k)][t] for k in j_to_k(j, net) ) == 0
                )
            end
        
        # unconnected nodes (should not exist so we warn)
        elseif isempty(i_to_j(j, net)) && isempty(j_to_k(j, net))
            @warn "Bus $j has no edges, setting loadbalcons to zero."
            m[:loadbalcons][j]["p"] = zeros(net.Ntimesteps)
            m[:loadbalcons][j]["q"] = zeros(net.Ntimesteps)

        # leaf nodes / sinks, flows in = draw out
        elseif !isempty(i_to_j(j, net)) && isempty(j_to_k(j, net))
            m[:loadbalcons][j]["p"] = @constraint(m, [t = 1:net.Ntimesteps],
                sum( pij[(i,j)][t] for i in i_to_j(j, net) )
                + pj[t] == 0
            )
            m[:loadbalcons][j]["q"] = @constraint(m, [t = 1:net.Ntimesteps],
                sum( qij[(i,j)][t] for i in i_to_j(j, net) )
                - shunt_susceptance * m[:vsqrd][j][t]
                + qj[t] + var == 0
            )
        
        # intermediate nodes
        else
            m[:loadbalcons][j]["p"] = @constraint(m, [t = 1:net.Ntimesteps],
                sum( pij[(i,j)][t] for i in i_to_j(j, net) )
                - sum( pij[(j,k)][t] for k in j_to_k(j, net) )
                + pj[t] == 0
            )
            m[:loadbalcons][j]["q"] = @constraint(m, [t = 1:net.Ntimesteps],
                sum( qij[(i,j)][t] for i in i_to_j(j, net) ) 
                - sum( qij[(j,k)][t] for k in j_to_k(j, net) ) 
                - shunt_susceptance * m[:vsqrd][j][t]
                + qj[t] + var == 0
            )
        end
    end

    nothing
end


function constrain_substation_voltage(m, net::Network{SinglePhase})
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


function constrain_KVL(m, net::Network{SinglePhase})
    w = m[:vsqrd]
    P = m[:pij]
    Q = m[:qij]
    l = m[:lij]
    m[:vcons] = Dict()
    for j in busses(net)
        for i in i_to_j(j, net)  # for radial network there is only one i in i_to_j
            if !( isa(net[(i,j)], CommonOPF.VoltageRegulator) )
                rᵢⱼ = rij_per_unit(i,j,net)
                xᵢⱼ = xij_per_unit(i,j,net)
                m[:vcons][j] = @constraint(m, [t = 1:net.Ntimesteps],
                    w[j][t] == w[i][t]
                        - 2*(rᵢⱼ * P[(i,j)][t] + xᵢⱼ * Q[(i,j)][t])
                        + (rᵢⱼ^2 + xᵢⱼ^2) * l[(i,j)][t]
                )
            else
                reg = net[(i,j)]
                if !ismissing(reg.vreg_pu)
                    m[:vcons][j] = @constraint(m, [t = 1:net.Ntimesteps],
                        w[j][t] == reg.vreg_pu^2
                    )
                else  # use turn ratio
                    m[:vcons][j] = @constraint(m, [t = 1:net.Ntimesteps],
                        w[j][t] == w[i][t] * reg.turn_ratio^2 
                    )
                end
            end
        end
    end
    nothing
end


function constrain_KVL_linear(m, net::Network{SinglePhase})
    w = m[:vsqrd]
    P = m[:pij]
    Q = m[:qij]
    m[:vcons] = Dict()
    for j in busses(net)
        for i in i_to_j(j, net)  # for radial network there is only one i in i_to_j
            if !( isa(net[(i,j)], CommonOPF.VoltageRegulator) )
                rᵢⱼ = rij_per_unit(i,j,net)
                xᵢⱼ = xij_per_unit(i,j,net)
                m[:vcons][j] = @constraint(m, [t = 1:net.Ntimesteps],
                    w[j][t] == w[i][t]
                        - 2*(rᵢⱼ * P[(i,j)][t] + xᵢⱼ * Q[(i,j)][t])
                )
            else
                reg = net[(i,j)]
                if !ismissing(reg.vreg_pu)
                    m[:vcons][j] = @constraint(m, [t = 1:net.Ntimesteps],
                        w[j][t] == reg.vreg_pu^2
                    )
                else  # use turn ratio
                    m[:vcons][j] = @constraint(m, [t = 1:net.Ntimesteps],
                        w[j][t] == w[i][t] * reg.turn_ratio^2 
                    )
                end
            end
        end
    end
    nothing
end


function constrain_cone(m, net::Network{SinglePhase})
    w = m[:vsqrd]
    P = m[:pij]
    Q = m[:qij]
    l = m[:lij]
    for j in busses(net)
        for i in i_to_j(j, net)  # for radial network there is only one i in i_to_j
            # # equivalent but maybe not handled as well in JuMP ?
            # @constraint(m, [t = 1:net.Ntimesteps],
            #     w[i][t] * l[(i,j)][t] ≥ P[(i,j)][t]^2 + Q[(i,j)][t]^2
            # )
            @constraint(m, [t = 1:net.Ntimesteps], 
                [w[i][t]/2, l[(i,j)][t], P[(i,j)][t], Q[(i,j)][t]] in JuMP.RotatedSecondOrderCone()
            )
        end
    end
end


function constrain_bilinear(m, net::Network{SinglePhase})
    w = m[:vsqrd]
    P = m[:pij]
    Q = m[:qij]
    l = m[:lij]
    for j in busses(net)
        for i in i_to_j(j, net)  # for radial network there is only one i in i_to_j
            @constraint(m, [t = 1:net.Ntimesteps],
                w[i][t] * l[(i,j)][t] == P[(i,j)][t]^2 + Q[(i,j)][t]^2
            )
        end
    end
end
