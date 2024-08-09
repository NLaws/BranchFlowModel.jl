@testset "NLP" begin
    dssfilepath = "data/ieee13/IEEE13Nodeckt_no_trfxs.dss"
    net = BranchFlowModel.CommonOPF.dss_to_Network(dssfilepath)

    net.v0 = 1.05
    net.Vbase = 4160 / sqrt(3)
    net.Sbase = 1e5
    net.Zbase = net.Vbase^2 / net.Sbase
    net.bounds.v_upper_mag = net.v0 * 1.1
    net.bounds.v_lower_mag = net.v0 * 0.8
    # net.bounds.s_upper_real =  1e8
    # net.bounds.s_lower_real = 0
    # net.bounds.i_upper_mag =  1e6
    net.bounds.i_lower_mag = 0


    m = Model(Ipopt.Optimizer)
    set_optimizer_attribute(m, "print_level", 0)

    build_model!(m, net, Unrelaxed)

    @objective(m, Min, 
        sum( 
            sum(real.(m[:i][t][i_j]) .* real.(m[:i][t][i_j])
            + imag.(m[:i][t][i_j]) .* imag.(m[:i][t][i_j])) 
            for t in 1:net.Ntimesteps, i_j in edges(net) 
        )
    )

    optimize!(m)
    
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

    vs = Dict(
        k => abs.(JuMP.value.(v))
        for (k,v) in m[:v][1]
    )

    i_ij = Dict(
        (i,j) => JuMP.value.(ii)
        for ((i,j), ii) in m[:i][1]
    )

    sj = Dict(
        k => JuMP.value.(sj)
        for (k,sj) in m[:Sj][1]
    )

    sij = Dict(
        (i,j) => JuMP.value.(ii)
        for ((i,j), ii) in m[:Sij][1]
    )

    # test power balance @ bus 671
    j = "671"
    net_power = sum( diag( 
        sij[(i,j)] - CPF.zij_per_unit(i,j,net) * i_ij[(i,j)] * BranchFlowModel.cj(i_ij[(i,j)])
    ) for i in CPF.i_to_j(j, net) ) - sum( 
        diag( sij[(j,k)] ) for k in CPF.j_to_k(j, net) )

    for phs in 1:3
        @test real(net_power[phs]) ≈ net[j, :kws, phs][1] * 1e3 / net.Sbase
        @test imag(net_power[phs]) ≈ net[j, :kvars, phs][1] * 1e3 / net.Sbase
    end


    # NOTE the OpenDSS model requires lots of iterations and setting load vminpu to 0.8 to ensure
    # that the OpenDSS model converges with constant power loads
    OpenDSS.Text.Command("Redirect $dssfilepath")

    OpenDSS.Solution.Solve()

    @test(OpenDSS.Solution.Converged() == true)
    @test(check_opendss_powers() == true)

    dss_voltages = dss_voltages_pu()

    errors = []
    for b in keys(vs)
        for (i, phsv) in enumerate(filter(v -> v != 0, vs[b]))
            @test abs(phsv - dss_voltages[b][i]) < 1e-3
            push!(errors, phsv - dss_voltages[b][i])
            # println("$b - $i  $(phsv - dss_voltages[b][i])")
        end
    end

end