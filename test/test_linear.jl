@testset "linear single phase model" begin

    dssfilepath = "data/singlephase38lines/master.dss"
    net = Network(dssfilepath)
    net.Sbase = 1e3  # put loads into kW and KVaR

    m = Model(HiGHS.Optimizer)
    build_bfm!(m, net, Linear)
    @objective(m, Min, 
        sum( 
            m[:p0][t] + m[:q0][t] for t in 1:net.Ntimesteps
        )
    )
    set_optimizer_attribute(m, "output_flag", false)

    optimize!(m)

    @test termination_status(m) == MOI.OPTIMAL

    # lossless model
    @test value.(m[:p0])[1] ≈ CPF.total_load_kw(net)[1] rtol=1e-7
    @test value.(m[:q0])[1] ≈ CPF.total_load_kvar(net)[1] rtol=1e-7
    
end


@testset "multiphase LinDistFlow Arnold 2016 validation" begin
    #=
    Validating against Arnold 2016 results, which requires:
    - matching the loads, capacitors, and lines in https://github.com/msankur/LinDist3Flow/tree/master/IEEE%20PES%202016/matlab/feeder13
    - changing the objective function in IEEE PES 2016/matlab/2016.01.05/V_Balance_Solver_Yapprox_20160105.m on lines 400 and 402 to Z = Z + (FV(:,:,k1)*X)'*(FV(:,:,k1)*X); and Z = Z + 0.5*(Fu*X)'*(Fu*X);, i.e. use the squared L2 norm rather than the L2 norm 
        - note that the optimal Qvar values reported in the paper used the L2 norm objective, but the addition of the sqrt to the NLobjective does not work with Ipopt, so I modified the matlab objective from msankure to match the paper objective and the objective herein, which results in different optimal Qvar values than those reported in the paper. This test compares the optimal Qvar values to those obtained from the Matlab model with the squared L2 norm objective.
    - removing the random additions to the a0 and a1 coefficients in IEEE PES 2016/matlab/2016.01.05/runSim_iter_1step.m (and running the same file)
    - increasing impedances by 1.25 (Section IV, paragraph 3)
    - use equ.s (3) & (4) of their paper for loads (the a0 and a1 terms)
    - allowing decisions for reactive power injection at busses 632, 675, 680, and 684
    - NOTE that the signs for loads are opposite to Arnold 2016 (they assume negative injections / positive loads)
    =#

    test_values = Dict(
        "632" => [-0.003589386844896, -0.012095150336674, 0.016058346153863],
        "680" => [0.001437239423406, 0.018476869715961, -0.019964857130616],
        "675" => [0.001294902379144, 0.019607280856338, -0.021111771515487],
        "684" => [0.001114839279833, 0.0, -0.019045362995131]
    )

    net = Network(joinpath("data", "ieee13", "IEEE13_Arnold_2016.dss"))
    net.Sbase = 5_000_000  # msankur has 5,000 kva * 1e3
    net.Vbase = 4160
    net.Zbase = net.Vbase^2 / net.Sbase
    net.bounds.v_lower_mag = 0.8
    net.bounds.v_upper_mag = 1.5

    for e in edges(net)
        net[e].rmatrix *= 1.25
        net[e].xmatrix *= 1.25
    end

    m = Model(Ipopt.Optimizer)
    set_optimizer_attribute(m, "print_level", 0)

    build_bfm!(m, net, Linear)

    a0 = 0.9
    a1 = 0.1

    pij = m[:pij]
    qij = m[:qij]

    # add loads with voltage dependency
    # similar constraints to constrain_linear_power_balance (started with copy/paste from there)
    for j in CPF.load_busses(net)
        sj = CPF.sj_per_unit(j, net)
        pj, qj = real(sj), imag(sj)

        for phs in CPF.phases_into_bus(net, j)

            delete(m, m[:power_balance_constraints][j][:real][phs])
            delete(m, m[:power_balance_constraints][j][:reactive][phs])

            ks_on_phs = [k for k in j_to_k(j, net) if phs in CPF.phases_into_bus(net, k)]

            if !isempty(ks_on_phs)  # mid node
                m[:power_balance_constraints][j][:real][phs] = @constraint(
                    m, [t in 1:net.Ntimesteps],
                    sum( pij[(i, j)][t][phs] for i in i_to_j(j, net) ) +
                    pj[phs][t] * (a0 + a1 * m[:vsqrd][j][t][phs]) - 
                    sum( pij[(j, k)][t][phs] for k in ks_on_phs ) == 0
                )
                m[:power_balance_constraints][j][:reactive][phs] = @constraint(
                    m, [t in 1:net.Ntimesteps],
                    sum( qij[(i, j)][t][phs] for i in i_to_j(j, net) ) +
                    qj[phs][t] * (a0 + a1 * m[:vsqrd][j][t][phs]) - 
                    sum( qij[(j, k)][t][phs] for k in ks_on_phs ) == 0
                )
            else  # leaf node
                m[:power_balance_constraints][j][:real][phs] = @constraint(
                    m, [t in 1:net.Ntimesteps],
                    sum( pij[(i, j)][t][phs] for i in i_to_j(j, net) ) + 
                    pj[phs][t] * (a0 + a1 * m[:vsqrd][j][t][phs]) == 0
                )
                m[:power_balance_constraints][j][:reactive][phs] = @constraint(
                    m, [t in 1:net.Ntimesteps],
                    sum( qij[(i, j)][t][phs] for i in i_to_j(j, net) ) + 
                    qj[phs][t] * (a0 + a1 * m[:vsqrd][j][t][phs]) == 0
                )
            end
        end
    end

    # add the var resources as variable var injections
    m[:Qvar] = Dict()
    Qresource_nodes = ["632", "675", "680", "684"]

    for j in Qresource_nodes

        m[:Qvar][j] = Dict()
        sj = CPF.sj_per_unit(j, net)
        pj, qj = real(sj), imag(sj)

        for phs in CPF.phases_into_bus(net, j)

            m[:Qvar][j][phs] = @variable(m, [1:net.Ntimesteps])
            delete(m, m[:power_balance_constraints][j][:reactive][phs])

            ks_on_phs = [k for k in j_to_k(j, net) if phs in CPF.phases_into_bus(net, k)]

            if !isempty(ks_on_phs)  # mid node 632 and 684 among the Qresource_nodes

                m[:power_balance_constraints][j][:reactive][phs] = @constraint(
                    m, [t in 1:net.Ntimesteps],
                    sum( qij[(i, j)][t][phs] for i in i_to_j(j, net) ) +
                    qj[phs][t] * (a0 + a1 * m[:vsqrd][j][t][phs]) - 
                    sum( qij[(j, k)][t][phs] for k in ks_on_phs ) +
                    m[:Qvar][j][phs][t] / net.Sbase == 0
                )

            else  # leaf node

                if j == "675"  # capacitor injecting 200 kvar on all three phases

                    m[:power_balance_constraints][j][:reactive][phs] = @constraint(
                        m, [t in 1:net.Ntimesteps],
                        sum( qij[(i, j)][t][phs] for i in i_to_j(j, net) ) + 
                        qj[phs][t] * (a0 + a1 * m[:vsqrd][j][t][phs]) +
                        200_000 / net.Sbase +
                        m[:Qvar][j][phs][t] / net.Sbase == 0
                    )

                else
                    m[:power_balance_constraints][j][:reactive][phs] = @constraint(
                        m, [t in 1:net.Ntimesteps],
                        sum( qij[(i, j)][t][phs] for i in i_to_j(j, net) ) + 
                        qj[phs][t] * (a0 + a1 * m[:vsqrd][j][t][phs]) +
                        m[:Qvar][j][phs][t] / net.Sbase == 0
                    )
                end
            end
        end
    end

    # capacitor on bus 611, only has phase 3
    j = "611"
    phs = 3
    delete(m, m[:power_balance_constraints][j][:reactive][phs])
    qj = imag(CPF.sj_per_unit(j, net))
    m[:power_balance_constraints][j][:reactive][phs] = 
        @constraint(m, [t in 1:net.Ntimesteps],
            sum( qij[(i, j)][t][phs] for i in i_to_j(j, net) ) +
            qj[phs][t] * (a0 + a1 * m[:vsqrd][j][t][phs]) +
            100_000 / net.Sbase == 0 
        )
    
    @objective(m, Min,
        0.5* sum( (m[:Qvar][b][phs][1] / net.Sbase)^2 for b in Qresource_nodes, phs in CPF.phases_into_bus(net, b)) +
        sum(
            (m[:vsqrd][b][1][phs1] - m[:vsqrd][b][1][phs2])^2
            for b in setdiff(busses(net), [net.substation_bus, "rg60"]), 
            (phs1, phs2) in [[1,2], [1,3], [2,3]] 
            if phs1 in CPF.phases_into_bus(net, b) && phs2 in CPF.phases_into_bus(net, b)
        )
    )

    optimize!(m)
    @test termination_status(m) == MOI.LOCALLY_SOLVED
        
    # the validation
    for b in Qresource_nodes, phs in CPF.phases_into_bus(net, b)
        @test -test_values[b][phs] ≈ value(m[:Qvar][b][phs][1])/net.Sbase rtol=1e-3
    end

end
