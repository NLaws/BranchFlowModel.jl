"""
    add_linear_variables(m, net::Network{MultiPhase})

Define real phase-vector variables for:
- `:vsqrd` bus voltage magnitude squared
- `:pij`, `:qij` branch power flows
- `:p`, `:q` for the slack bus net power injection
"""
function add_linear_variables(m, net::Network{MultiPhase})
    net.var_names = [:vsqrd, :pij, :qij, :p, :q]

    m[:vsqrd] = multiphase_bus_variable_container()
    m[:p] = multiphase_bus_variable_container()
    m[:q] = multiphase_bus_variable_container()
    m[:pij] = multiphase_edge_variable_container()
    m[:qij] = multiphase_edge_variable_container()

    for b in busses(net), t in 1:net.Ntimesteps
        if b == net.substation_bus

            # TODO allow for time-varying source voltage
            m[:vsqrd][b][t] = abs.(substation_voltage_squared(net))

            m[:p][b][t] =  @variable(m, [1:3], base_name="p_$(b)_$(t)")
            m[:q][b][t] =  @variable(m, [1:3], base_name="q_$(b)_$(t)")

        else
            m[:vsqrd][b][t] = @variable(m, 
                [phs in phases_into_bus(net, b)], 
                lower_bound=0, 
                base_name="vsqrd_$(b)_$(t)"
            )
        end
    end

    for e in edges(net), t in 1:net.Ntimesteps
        m[:pij][e][t] = @variable(m, 
            [phs in phases_into_bus(net, e[2])], 
            base_name="pij_$(e)_$(t)"
        )
        m[:qij][e][t] = @variable(m, 
            [phs in phases_into_bus(net, e[2])], 
            base_name="qij_$(e)_$(t)"
        )
    end

    net.var_info[:vsqrd] = CommonOPF.VariableInfo(
        :vsqrd,
        "voltage magnitude squared",
        CommonOPF.VoltUnit,
        (CommonOPF.BusDimension, CommonOPF.TimeDimension, CommonOPF.PhaseDimension)
    )

    net.var_info[:pj] = CommonOPF.VariableInfo(
        :pj,
        "net bus injection real power on bus j",
        CommonOPF.RealPowerUnit,
        (CommonOPF.BusDimension, CommonOPF.TimeDimension, CommonOPF.PhaseDimension)
    )

    net.var_info[:qj] = CommonOPF.VariableInfo(
        :qj,
        "net bus injection reactive power on bus j",
        CommonOPF.ReactivePowerUnit,
        (CommonOPF.BusDimension, CommonOPF.TimeDimension, CommonOPF.PhaseDimension)
    )

    net.var_info[:pij] = CommonOPF.VariableInfo(
        :pij,
        "sending end real power from bus i to j",
        CommonOPF.RealPowerUnit,
        (CommonOPF.EdgeDimension, CommonOPF.TimeDimension, CommonOPF.PhaseDimension)
    )

    net.var_info[:qij] = CommonOPF.VariableInfo(
        :qij,
        "sending end reactive power from bus i to j",
        CommonOPF.ReactivePowerUnit,
        (CommonOPF.EdgeDimension, CommonOPF.TimeDimension, CommonOPF.PhaseDimension)
    )

    # TODO net.bounds

    nothing
end


"""
    constrain_linear_power_balance(m, net::Network{MultiPhase})
``
p_{ij,\\phi} + p_{j,\\phi} = \\sum{k:j\\rightarrow k} p_{jk,\\phi} 
\\quad \\forall j \\in \\mathcal{N},
\\forall \\phi \\in [1,2,3]
``

``
q_{ij,\\phi} + q_{j,\\phi} = \\sum{k:j\\rightarrow k} q_{jk,\\phi} 
\\quad \\forall j \\in \\mathcal{N}, 
\\forall \\phi \\in [1,2,3]
``
"""
function constrain_linear_power_balance(m, net::Network{MultiPhase})
    p0 = m[:p]
    q0 = m[:q]
    pij = m[:pij]
    qij = m[:qij]

    for j in busses(net)
        pj, qj = sj_per_unit(j, net)

        if j == net.substation_bus   # include the slack power variables

            for phs in [1,2,3]  # TODO can vectorize constraints across phs?
                ks_on_phs = [k for k in j_to_k(j, net) if phs in phases_into_bus(net, k)]
                @constraint(m, [t in 1:net.Ntimesteps],
                    pj[phs][t] + p0[j][t][phs] - sum( pij[(j, k)][t][phs] for k in ks_on_phs ) == 0
                )
                @constraint(m, [t in 1:net.Ntimesteps],
                    qj[phs][t] + q0[j][t][phs] - sum( qij[(j, k)][t][phs] for k in ks_on_phs ) == 0
                )
            end     

        elseif isempty(i_to_j(j, net)) && isempty(j_to_k(j, net))  # unconnected nodes
            @warn "Bus $j has no edges in or out; constraining pj and qj to zero."
            @constraint(m, [phs in 1:3, t in 1:net.Ntimesteps],
                pj[phs][t] == 0
            )
            @constraint(m, [phs in 1:3, t in 1:net.Ntimesteps],
                qj[phs][t] == 0
            )

        else  # mid and leaf nodes

            for phs in phases_into_bus(net, j)
                ks_on_phs = [k for k in j_to_k(j, net) if phs in phases_into_bus(net, k)]
                if !isempty(ks_on_phs)  # mid node
                    @constraint(m, [t in 1:net.Ntimesteps],
                        sum( pij[(i, j)][t][phs] for i in i_to_j(j, net) ) +
                        pj[phs][t] - sum( pij[(j, k)][t][phs] for k in ks_on_phs ) == 0
                    )
                    @constraint(m, [t in 1:net.Ntimesteps],
                        sum( qij[(i, j)][t][phs] for i in i_to_j(j, net) ) +
                        qj[phs][t] - sum( qij[(j, k)][t][phs] for k in ks_on_phs ) == 0
                    )
                else  # leaf node
                    @constraint(m, [t in 1:net.Ntimesteps],
                        sum( pij[(i, j)][t][phs] for i in i_to_j(j, net) ) + pj[phs][t] == 0
                    )
                    @constraint(m, [t in 1:net.Ntimesteps],
                        sum( qij[(i, j)][t][phs] for i in i_to_j(j, net) ) + qj[phs][t] == 0
                    )
                end
            end

        end
    end
    nothing
end


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


"""
    constrain_KVL_linear(m, net::Network{MultiPhase})

``
    \\boldsymbol{w}_j = \\boldsymbol{w}_i + \\boldsymbol{M}_{P,ij} \\boldsymbol{P}_{ij} 
    + \\boldsymbol{M}_{Q,ij} \\boldsymbol{Q}_{ij}
``

See [`MPij`](@ref) and [`MQij`](@ref) for their definitions.
"""
function constrain_KVL_linear(m, net::Network{MultiPhase})
    w = m[:vsqrd]
    p = m[:pij]
    q = m[:qij]
    m[:kvl_constraints] = Dict()

    for j in busses(net)  # substation_bus in here but has empty i_to_j(j, net)
        for i in i_to_j(j, net)  # for radial network there is only one i in i_to_j
            i_j = (i, j)

            if !( isa(net[(i,j)], CommonOPF.VoltageRegulator) )
                
                MP = MPij(i,j,net)
                MQ = MQij(i,j,net)
                m[:kvl_constraints][i_j] = Dict()
                for phs in phases_into_bus(net, j)
                    m[:kvl_constraints][i_j][phs] = @constraint(m, [t in 1:net.Ntimesteps],
                        w[j][t][phs] == w[i][t][phs]
                            + sum(MP[phs,k] * p[i_j][t][k] for k=phases_into_bus(net, j))
                            + sum(MQ[phs,k] * q[i_j][t][k] for k=phases_into_bus(net, j))
                    )
                end

            else
                reg = net[(i,j)]
                if !ismissing(reg.vreg_pu)
                    m[:kvl_constraints][i_j] = @constraint(m, [t in 1:net.Ntimesteps],
                        w[j][t] .== reg.vreg_pu.^2
                    )
                else  # use turn ratio
                    m[:kvl_constraints][i_j] = @constraint(m, [t in 1:net.Ntimesteps],
                        w[j][t] .== w[i][t] * reg.turn_ratio^2 
                    )
                end
            end
        end
    end

    # document the constraints
    e = edges(net)[1]
    phs = collect(keys(m[:kvl_constraints][e][1]))[1]
    c = m[:kvl_constraints][e][1][phs]  # time step 1
    net.constraint_info[:kvl_constraints] = CommonOPF.ConstraintInfo(
        :kvl_constraints,
        "Kirchoff's Voltage Law linear version",
        MOI.get(m, MOI.ConstraintSet(), c),
        (CommonOPF.EdgeDimension, CommonOPF.TimeDimension, CommonOPF.PhaseDimension),
    )
    nothing
end
