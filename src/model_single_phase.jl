"""
    build_model!(m::JuMP.AbstractModel, p::Inputs)

Add variables and constraints to `m` using the values in `p`. Calls the following functions:
```julia
add_variables(m, p)
constrain_power_balance(m, p)
constrain_substation_voltage(m, p)
constrain_KVL(m, p)
constrain_cone(m, p)
constrain_loads(m, p)
```
"""
function build_model!(m::JuMP.AbstractModel, p::Inputs)

    add_variables(m, p)
    constrain_power_balance(m, p)
    constrain_substation_voltage(m, p)
    constrain_KVL(m, p)
    constrain_loads(m, p)
    if p.relaxed
        constrain_cone(m, p)
    else
        constrain_psd(m, p)
    end
end


function add_variables(m, p::Inputs)
    T = 1:p.Ntimesteps
    # bus injections
    @variables m begin
        p.P_lo_bound <= Pj[p.busses, T] <= p.P_up_bound
        p.Q_lo_bound <= Qj[p.busses, T] <= p.Q_up_bound
    end

    # voltage squared
    @variable(m, p.v_lolim^2 <= vsqrd[p.busses, T] <= p.v_uplim^2 ) 
    
    # line flows, power sent from i to j
    @variable(m, p.P_lo_bound <= Pij[p.edge_keys, T] <= p.P_up_bound )
    
    @variable(m, p.Q_lo_bound <= Qij[p.edge_keys, T] <= p.Q_up_bound )

    # current squared (non-negative)
    @variable(m, 0 <= lij[p.edge_keys, T])

    # TODO better indexing for lines and attributes
    @constraint(m, [((i,j), edge) in zip(p.edges, p.edge_keys), t in T],
        lij[edge, t] <= p.Isquared_up_bounds[get_ijlinecode(i,j,p)]
    )

    nothing
end


"""
    function constrain_power_balance(m, p::Inputs)

Pij in - losses == sum of line flows out + net injection
NOTE: using sum over Pij for future expansion to mesh grids
i -> j -> k
"""
function constrain_power_balance(m, p::Inputs)
    Pj = m[:Pj]
    Qj = m[:Qj]
    Pij = m[:Pij]
    Qij = m[:Qij]
    lij = m[:lij]
    m[:loadbalcons] = Dict()
    # TODO change Pj and Qj to expressions, make P₀ and Q₀ dv's, which will reduce # of variables
    # by (Nnodes - 1)*8760 and number of constraints by 6*(Nnodes - 1)*8760
    for j in p.busses
        if isempty(i_to_j(j, p)) && !isempty(j_to_k(j, p)) # source nodes, injection = flows out
            substation_Pload = zeros(p.Ntimesteps)
            if j in keys(p.Pload)
                substation_Pload = p.Pload[j] / p.Sbase
            end
            pcon = @constraint(m,  [t in 1:p.Ntimesteps],
                Pj[j,t] - substation_Pload[t] - sum( Pij[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
            substation_Qload = zeros(p.Ntimesteps)
            if j in keys(p.Qload)
                substation_Qload = p.Qload[j] / p.Sbase
            end
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                Qj[j,t] - substation_Qload[t] - sum( Qij[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
        elseif isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # unconnected nodes
            @warn "Bus $j has no edges, setting Pj and Qj to zero."
            pcon = @constraint(m, [t in 1:p.Ntimesteps],
                Pj[j,t] == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                Qj[j,t] == 0
            )
        elseif !isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # leaf nodes / sinks, flows in = draw out
            pcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Pij[string(i*"-"*j), t] for i in i_to_j(j, p) ) -
                sum( lij[string(i*"-"*j), t] * rij(i,j,p) for i in i_to_j(j, p) ) 
                + Pj[j, t] == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Qij[string(i*"-"*j), t] for i in i_to_j(j, p) ) -
                sum( lij[string(i*"-"*j), t] * xij(i,j,p) for i in i_to_j(j, p) ) + 
                Qj[j, t] == 0
            )
        else
            pcon =  @constraint(m, [t in 1:p.Ntimesteps],
                sum( Pij[string(i*"-"*j), t] for i in i_to_j(j, p) ) - 
                sum( lij[string(i*"-"*j), t] * rij(i,j,p) for i in i_to_j(j, p) ) +
                Pj[j,t] - sum( Pij[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Qij[string(i*"-"*j), t] for i in i_to_j(j, p) ) 
                - sum( lij[string(i*"-"*j), t] * xij(i,j,p) for i in i_to_j(j, p) ) +
                Qj[j,t] - sum( Qij[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
        end
        m[:loadbalcons][j] = Dict("p" => pcon, "q" => qcon)
    end
    # TODO Farivar and Low have b*v and q*v in these equations for shunts? neglected for now

    nothing
end


function constrain_substation_voltage(m, p::Inputs)
    if typeof(p.v0) <: Real
        @constraint(m, con_substationV[t in 1:p.Ntimesteps],
            m[:vsqrd][p.substation_bus, t] == p.v0^2
        )
    else  # vector of time
        @constraint(m, con_substationV[t in 1:p.Ntimesteps],
            m[:vsqrd][p.substation_bus, t] == p.v0[t]^2
        )
    end
    nothing
end


function constrain_KVL(m, p::Inputs)
    w = m[:vsqrd]
    P = m[:Pij]
    Q = m[:Qij]
    l = m[:lij]
    m[:vcons] = Dict()
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            if !( (i,j) in keys(p.regulators) )
                i_j = string(i*"-"*j)
                rᵢⱼ = rij(i,j,p)
                xᵢⱼ = xij(i,j,p)
                m[:vcons][j] = @constraint(m, [t in 1:p.Ntimesteps],
                    w[j,t] == w[i,t]
                        - 2*(rᵢⱼ * P[i_j,t] + xᵢⱼ * Q[i_j,t])
                        + (rᵢⱼ^2 + xᵢⱼ^2) * l[i_j, t]
                )
            else
                if has_vreg(p, j)
                    m[:vcons][j] = @constraint(m, [t in 1:p.Ntimesteps],
                        w[j,t] == p.regulators[(i,j)][:vreg]^2
                    )
                else  # default turn_ratio is 1.0
                    m[:vcons][j] = @constraint(m, [t in 1:p.Ntimesteps],
                        w[j,t] == w[i,t] * p.regulators[(i,j)][:turn_ratio]^2 
                    )
                end
            end
        end
    end
    nothing
end


function constrain_cone(m, p::Inputs)
    w = m[:vsqrd]
    P = m[:Pij]
    Q = m[:Qij]
    l = m[:lij]
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            i_j = string(i*"-"*j)
            # # equivalent but maybe not handled as well in JuMP ?
            # @constraint(m, [t in 1:p.Ntimesteps],
            #     w[i,t] * l[i_j, t] ≥ P[i_j,t]^2 + Q[i_j,t]^2
            # )
            @constraint(m, [t in 1:p.Ntimesteps], 
                [w[i,t]/2, l[i_j, t], P[i_j,t], Q[i_j,t]] in JuMP.RotatedSecondOrderCone()
            )
        end
    end
end


function constrain_psd(m, p::Inputs)
    w = m[:vsqrd]
    P = m[:Pij]
    Q = m[:Qij]
    l = m[:lij]
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            i_j = string(i*"-"*j)
            # need to use PSD?
            @constraint(m, [t in 1:p.Ntimesteps],
                w[i,t] * l[i_j, t] == P[i_j,t]^2 + Q[i_j,t]^2
            )
        end
    end
end


"""
    constrain_loads(m, p::Inputs)

- set net injections Pj/Qj to negative of Inputs.Pload/Qload, which are normalized by Sbase when creating Inputs
- keys of P/Qload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus
"""
function constrain_loads(m, p::Inputs)
    Pj = m[:Pj]
    Qj = m[:Qj]
    m[:injectioncons] = Dict()
    for j in setdiff(p.busses, [p.substation_bus])
        m[:injectioncons][j] = Dict()
        if j in keys(p.Pload)
            con = @constraint(m, [t in 1:p.Ntimesteps],
                Pj[j,t] == -p.Pload[j][t] / p.Sbase
            )
        else
            con = @constraint(m, [t in 1:p.Ntimesteps],
                Pj[j,t] == 0
            )
        end
        m[:injectioncons][j]["p"] = con
        if j in keys(p.Qload)
            con = @constraint(m, [t in 1:p.Ntimesteps],
                Qj[j,t] == -p.Qload[j][t] / p.Sbase
            )
        else
            con = @constraint(m, [t in 1:p.Ntimesteps],
                Qj[j,t] == 0
            )
        end
        m[:injectioncons][j]["q"] = con
    end
    nothing
end
