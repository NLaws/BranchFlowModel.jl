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