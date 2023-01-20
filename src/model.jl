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
    constrain_cone(m, p)
    constrain_loads(m, p)

end


function add_variables(m, p::Inputs)
    T = 1:p.Ntimesteps
    # bus injections
    @variables m begin
        p.P_lo_bound <= Pⱼ[p.busses, T] <= p.P_up_bound
        p.Q_lo_bound <= Qⱼ[p.busses, T] <= p.Q_up_bound
    end

    # voltage squared
    @variable(m, p.v_lolim^2 <= vsqrd[p.busses, T] <= p.v_uplim^2 ) 
    
    ij_edges = [string(i*"-"*j) for j in p.busses for i in i_to_j(j, p)]

    # line flows, power sent from i to j
    @variable(m, p.P_lo_bound <= Pᵢⱼ[ij_edges, T] <= p.P_up_bound )
    
    @variable(m, p.Q_lo_bound <= Qᵢⱼ[ij_edges, T] <= p.Q_up_bound )

    # current squared (non-negative)
    @variable(m, 0 <= lᵢⱼ[ij_edges, T])

    # TODO better indexing for lines and attributes
    @constraint(m, [((i,j), edge) in zip(p.edges, ij_edges), t in T],
        lᵢⱼ[edge, t] <= p.Isqaured_up_bounds[get_ijlinecode(i,j,p)]
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
    Pⱼ = m[:Pⱼ]
    Qⱼ = m[:Qⱼ]
    Pᵢⱼ = m[:Pᵢⱼ]
    Qᵢⱼ = m[:Qᵢⱼ]
    lᵢⱼ = m[:lᵢⱼ]
    m[:loadbalcons] = Dict()
    # TODO change Pⱼ and Qⱼ to expressions, make P₀ and Q₀ dv's, which will reduce # of variables
    # by (Nnodes - 1)*8760 and number of constraints by 6*(Nnodes - 1)*8760
    for j in p.busses
        if isempty(i_to_j(j, p)) && !isempty(j_to_k(j, p)) # source nodes, injection = flows out
            pcon = @constraint(m,  [t in 1:p.Ntimesteps],
                Pⱼ[j,t] - sum( Pᵢⱼ[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] - sum( Qᵢⱼ[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
        elseif isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # unconnected nodes
            @warn "Bus $j has no edges, setting Pⱼ and Qⱼ to zero."
            pcon = @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,t] == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == 0
            )
        elseif !isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # leaf nodes / sinks, flows in = draw out
            pcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Pᵢⱼ[string(i*"-"*j), t] for i in i_to_j(j, p) ) -
                sum( lᵢⱼ[string(i*"-"*j), t] * rij(i,j,p) for i in i_to_j(j, p) ) 
                + Pⱼ[j, t] == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Qᵢⱼ[string(i*"-"*j), t] for i in i_to_j(j, p) ) -
                sum( lᵢⱼ[string(i*"-"*j), t] * xij(i,j,p) for i in i_to_j(j, p) ) + 
                Qⱼ[j, t] == 0
            )
        else
            pcon =  @constraint(m, [t in 1:p.Ntimesteps],
                sum( Pᵢⱼ[string(i*"-"*j), t] for i in i_to_j(j, p) ) - 
                sum( lᵢⱼ[string(i*"-"*j), t] * rij(i,j,p) for i in i_to_j(j, p) ) +
                Pⱼ[j,t] - sum( Pᵢⱼ[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Qᵢⱼ[string(i*"-"*j), t] for i in i_to_j(j, p) ) 
                - sum( lᵢⱼ[string(i*"-"*j), t] * xij(i,j,p) for i in i_to_j(j, p) ) +
                Qⱼ[j,t] - sum( Qᵢⱼ[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
        end
        m[:loadbalcons][j] = Dict("p" => pcon, "q" => qcon)
    end
    # TODO Farivar and Low have b*v and q*v in these equations for shunts? neglected for now

    nothing
end


function constrain_substation_voltage(m, p::Inputs)
    # @info "constrain_substation_voltage"
    @constraint(m, con_substationV[t in 1:p.Ntimesteps],
       m[:vsqrd][p.substation_bus, t] == p.v0^2
    )
    nothing
end


function constrain_KVL(m, p::Inputs)
    w = m[:vsqrd]
    P = m[:Pᵢⱼ]
    Q = m[:Qᵢⱼ]
    l = m[:lᵢⱼ]
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            i_j = string(i*"-"*j)
            rᵢⱼ = rij(i,j,p)
            xᵢⱼ = xij(i,j,p)
            vcon = @constraint(m, [t in 1:p.Ntimesteps],
                w[j,t] == w[i,t]
                    - 2*(rᵢⱼ * P[i_j,t] + xᵢⱼ * Q[i_j,t])
                    + (rᵢⱼ^2 + xᵢⱼ^2) * l[i_j, t]
            )
        end
    end
    nothing
end


function constrain_cone(m, p::Inputs)
    w = m[:vsqrd]
    P = m[:Pᵢⱼ]
    Q = m[:Qᵢⱼ]
    l = m[:lᵢⱼ]
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


"""
    constrain_loads(m, p::Inputs)

- set loads to negative of Inputs.Pload, which are normalized by Sbase when creating Inputs
- keys of Pload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus
"""
function constrain_loads(m, p::Inputs)
    Pⱼ = m[:Pⱼ]
    Qⱼ = m[:Qⱼ]
    
    for j in p.busses
        if j in keys(p.Pload)
            @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,t] == -p.Pload[j][t] / p.Sbase
            )
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,t] == 0
            )
        end
        if j in keys(p.Qload)
            @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == -p.Qload[j][t] / p.Sbase
            )
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == 0
            )
        end
    end
    nothing
end
