using BranchFlowModel
using Test
using Random
using ECOS
using JuMP
using OpenDSSDirect
using Ipopt

Random.seed!(42)


@testset "BranchFlowModel.jl" begin


@testset "input construction" begin
    # then try test Inputs
    # IEEE 4 bus model for first tests?
    # need to test convergence using this module in UL with LL DER

    # lines
    edges = [("0", "1"), ("1", "2")]
    linecodes = ["0_1", "1_2"]
    linelengths = [100.0, 200.0]
    phases = repeat([[1]], length(edges))
    Zdict = Dict{String, Dict{String, Any}}(
        lc => Dict("rmatrix" => 0.001, "xmatrix" => 0.001)
        for lc in linecodes
    )

    # busses and loads
    busses = ["0", "1", "2"]
    substation_bus = "0"
    load_busses = setdiff(busses, [substation_bus])
    Pload = Dict(b => [10] for b in load_busses)
    Qload = Dict(b => [1]  for b in load_busses)

    # network
    Sbase = 1
    Vbase = 1

    Ntimesteps = 1

    # p = Inputs(
    #     edges,
    #     linecodes,
    #     linelengths,
    #     busses,
    #     phases,
    #     substation_bus,
    #     Pload,
    #     Qload,
    #     Sbase,
    #     Vbase,
    #     Ibase,
    #     Zdict,
    #     v0,
    #     v_lolim, 
    #     v_uplim,
    #     Zbase,
    #     Ntimesteps,
    #     0.1,  # power factor
    #     length(busses),  # Nnodes
    # )


    d = open(joinpath("data", "singlephase38lines", "master.dss")) do io  # 
        parse_dss(io)  # method from PowerModelsDistribution
    end
end

# using data taken from Andrianesis, Caramanis LMV paper 2019
# TODO don't use random loads
@testset "single phase 38-nodes 3 time steps" begin
    T = 3
    loadnodes = ["3", "5", "36", "9", "10", "11", "12", "13", "15", "17", "18", "19", "22", "25", 
                "27", "28", "30", "31", "32", "33", "34", "35"]

    loads = rand(length(loadnodes), T) * 1e2

    Pload = Dict(k =>     loads[indexin([k], loadnodes)[1], :] for k in loadnodes)
    Qload = Dict(k => 0.1*loads[indexin([k], loadnodes)[1], :] for k in loadnodes)

    Sbase = 1e6
    Vbase = 12.5e3

    inputs = Inputs(
        joinpath("data", "singlephase38lines", "master.dss"), 
        "0";
        Pload=Pload, 
        Qload=Qload,
        Sbase=Sbase, 
        Vbase=Vbase, 
        v0 = 1.00,
        v_uplim = 1.05,
        v_lolim = 0.95,
        Ntimesteps = T
    );

    m = Model(ECOS.Optimizer)
    build_model!(m, inputs)

    set_optimizer_attribute(m, "maxit", 10000)
    set_optimizer_attribute(m, "verbose", 0)
    #= ECOS OPTIONS
    gamma          # scaling the final step length
    delta          # regularization parameter
    eps            # regularization threshold
    feastol        # primal/dual infeasibility tolerance
    abstol         # absolute tolerance on duality gap
    reltol         # relative tolerance on duality gap
    feastol_inacc  # primal/dual infeasibility relaxed tolerance
    abstol_inacc   # absolute relaxed tolerance on duality gap
    reltol_inacc   # relative relaxed tolerance on duality gap
    nitref         # number of iterative refinement steps
    maxit          # maximum number of iterations
    verbose        # verbosity bool for PRINTLEVEL < 3
    =#
    # can add objective here
    @objective(m, Min, sum(m[:Pⱼ]["0", t] for t=1:T))
    optimize!(m)

    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]

end


@testset "ieee13 unbalanced MultiPhase" begin
    T = 3

    # make the dss solution to compare
    dss("Redirect data/ieee13/IEEE13Nodeckt.dss")
    @test(OpenDSSDirect.Solution.Converged() == true)

    dss_voltages = Dict(
        k => v for (k,v) in zip(
            OpenDSSDirect.Circuit.AllBusNames(), OpenDSSDirect.Circuit.AllBusVMag()
        )
    )

    inputs = Inputs(
        joinpath("data", "ieee13", "IEEE13Nodeckt.dss"), 
        "0";
        Sbase=5_000_000, 
        Vbase=4160, 
        v0 = 1.00,
        v_uplim = 1.05,
        v_lolim = 0.95,
        Ntimesteps = T
    );

    m = Model(Ipopt.Optimizer)


end

end  # all tests
