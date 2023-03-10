using BranchFlowModel
using Test
using Random
using ECOS
using JuMP
using OpenDSSDirect
using SCS
using LinearAlgebra
using COSMO
using CSDP
using Ipopt

# # hack for local testing
# using Pkg
# Pkg.activate("..")
# using BranchFlowModel
# Pkg.activate(".")

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
    @objective(m, Min, sum(m[:P???]["0", t] for t=1:T))
    optimize!(m)

    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]

end


function dss_voltages_pu()
    d = Dict()
    for b in OpenDSSDirect.Circuit.AllBusNames() 
        OpenDSSDirect.Circuit.SetActiveBus(b)
        d[b] = OpenDSSDirect.Bus.puVmagAngle()[1:2:end]
    end
    return d
end


@testset "ieee13 unbalanced MultiPhase" begin

    # make the dss solution to compare
    dss("Redirect data/ieee13/IEEE13Nodeckt.dss")
    @test(OpenDSSDirect.Solution.Converged() == true)

    dss_voltages = dss_voltages_pu()

    p = Inputs(
        joinpath("data", "ieee13", "IEEE13Nodeckt.dss"), 
        "rg60";
        Sbase=1_000_000, 
        Vbase=4160, 
        v0 = dss_voltages["rg60"],  # openDSS rg60 values
        v_uplim = 1.06,
        v_lolim = 0.94,
        Ntimesteps = 1
    );
    # p.Isqaured_up_bounds = Dict(
    #     lc => 100 for lc in Set(p.linecodes)
    # )
    p.P_lo_bound = -10
    p.Q_lo_bound = -10
    p.P_up_bound = 10
    p.Q_up_bound = 10

    m = Model(CSDP.Optimizer)

    # m = Model(COSMO.Optimizer)

    # m = Model(SCS.Optimizer)

    build_model!(m,p)

    ij_edges = [string(i*"-"*j) for j in p.busses for i in i_to_j(j, p)];

    @objective(m, Min, 
        sum( sum(real.(diag(m[:l][t][i_j]))) for t in 1:p.Ntimesteps, i_j in  ij_edges)
    )

    optimize!(m)
    
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]

    @test_nowarn(check_rank_one(m,p))

    vs = Dict(k => real.(diag(v)) for (k,v) in voltage_values_by_time_bus(m,p)[1])

    # vs["632"]
    # dss_voltages["632"]

    # TODO why are BFM voltages so much higher than openDSS ?

    #=
        vs = Dict(k => real.(diag(v)) for (k,v) in voltage_values_by_time_bus(m,p)[1])
        Dict{String, Vector{Float64}} with 14 entries:
        "671"  => [1.02746, 1.02239, 1.0332]
        "680"  => [1.02746, 1.02239, 1.0332]
        "652"  => [1.02421, 0.0, 0.0]
        "634"  => [1.04153, 1.02472, 1.04205]
        "675"  => [1.02545, 1.02258, 1.03169]
        "rg60" => [1.05603, 1.03739, 1.05605]
        "611"  => [0.0, 0.0, 1.031]
        "645"  => [0.0, 1.02191, 1.04052]
        "632"  => [1.04258, 1.02554, 1.04262]
        "633"  => [1.04154, 1.02473, 1.04206]
        "684"  => [1.02613, 0.0, 1.03208]
        "692"  => [1.02746, 1.02239, 1.03319]
        "670"  => [1.0374, 1.02425, 1.03881]
        "646"  => [0.0, 1.02133, 1.03973]

        dss_voltages
        Dict{Any, Any} with 16 entries:
        "671"       => [0.982795, 1.04028, 0.964876]
        "680"       => [0.982795, 1.04028, 0.964876]
        "634"       => [0.987157, 1.00842, 0.982446]
        "652"       => [0.975332]
        "675"       => [0.97627, 1.04263, 0.962942]
        "650"       => [0.999911, 0.999971, 0.999931]
        "rg60"      => [1.05603, 1.03739, 1.05605]
        "611"       => [0.96083]
        "645"       => [1.01973, 1.00227]
        "632"       => [1.01434, 1.02895, 1.00419]
        "633"       => [1.0113, 1.02702, 1.00154]
        "684"       => [0.980872, 0.962846]
        "sourcebus" => [0.999974, 0.999994, 0.99995]
        "692"       => [0.982795, 1.04028, 0.964876]
        "670"       => [1.00402, 1.03187, 0.989743]
        "646"       => [1.01801, 1.00024]

        OpenDSSDirect.Circuit.SetActiveElement("Line.650632")

        OpenDSSDirect.CktElement.Powers()
            6-element Vector{ComplexF64}:
            1251.8855692490793 + 684.8666548773817im
            972.7880484340428 + 379.0983497071919im
            1342.1833785467552 + 672.0212464201516im
            -1230.1864192003165 - 604.4012027674321im
            -975.8849719090533 - 346.4379601073651im
            -1300.0498332554164 - 589.1374409448683im


        p.Sbase * value.(m[:Sj][1][p.substation_bus])/1e3
            3-element Vector{ComplexF64}:
            1282.0316166722455 + 736.6559894012869im
            933.7883751270138 + 613.6212174838729im
            1288.256272184194 + 865.9247220594767im

        OpenDSSDirect.CktElement.Currents()
            6-element Vector{ComplexF64}:
            493.51471731379274 - 270.13263500932226im
            -327.041545966684 - 261.97567100768674im
            -34.96614764181868 + 590.759580432783im
            -493.5146926040484 + 270.1338402015706im
            327.04257828741856 + 261.9750513219735im
            34.96509469955981 - 590.7601622800512im

        p.Ibase * sqrt.(value.(m[:l][1][p.substation_bus * "-632"]))
            3??3 Matrix{ComplexF64}:
            194.322-0.0im      89.5948-144.986im  93.4316+175.816im
            89.5948+144.986im  149.485-0.0im      88.1002-150.773im
            93.4316-175.816im  88.1002+150.773im  203.995-0.0im

        Z=(OpenDSSDirect.CktElement.YPrim())^-1

        julia> Z[1:3,1:3]
        3??3 Matrix{ComplexF64}:
        0.030064-2.83202e6im  0.0145558-7.7237e5im   0.0144071-7.7237e5im
        0.0137765-7.7237e5im    0.032481-2.83202e6im  0.0142786-7.7237e5im
        0.0135392-7.7237e5im   0.0148967-7.7237e5im   0.0323955-2.83202e6im

        julia> zij("rg60","632",p)
        3??3 Matrix{ComplexF64}:
        0.0379213+0.1114im     0.0170728+0.0549065im  0.0172917+0.0463591im
        0.0170728+0.0549065im  0.0369363+0.114672im   0.0167992+0.0421238im
        0.0172917+0.0463591im  0.0167992+0.0421238im  0.0373631+0.113249im

        similar R values, but openDSS has huge X values ?

        OpenDSSDirect.Lines.RMatrix() * OpenDSSDirect.Lines.Length() ./ real.(zij("rg60","632",p))
            3??3 Matrix{Float64}:
            3.46112  3.46112  3.46112
            3.46112  3.46112  3.46112
            3.46112  3.46112  3.46112

        OpenDSSDirect.Lines.XMatrix() * OpenDSSDirect.Lines.Length() ./ imag.(zij("rg60","632",p))
            3??3 Matrix{Float64}:
            3.46112  3.46112  3.46112
            3.46112  3.46112  3.46112
            3.46112  3.46112  3.46112

        add caps in manually

        # p.Qload["611"][3][1] -= 100_000  # eyeballing BFM vs DSS voltages I turned this off
        p.Qload["675"][1][1] -= 200_000
        p.Qload["675"][2][1] -= 200_000
        p.Qload["675"][3][1] -= 200_000
        
        what is up with openDSS X values so high?
    =#
end

@testset "ieee13 balanced SinglePhase" begin

    # make the dss solution to compare
    dss("Redirect data/ieee13_makePosSeq/Master.dss")
    dss("Solve")
    @test(OpenDSSDirect.Solution.Converged() == true)

    dss_voltages = dss_voltages_pu()

    vbase = 4160/sqrt(3)
    p = Inputs(
        joinpath("data", "ieee13_makePosSeq", "Master.dss"), 
        "rg60";
        Sbase=1_000_000, 
        Vbase=vbase, 
        v0 = dss_voltages["rg60"][1],  # openDSS rg60 values
        v_uplim = 1.05,
        v_lolim = 0.95,
        Ntimesteps = 1
    );
    p.relaxed = false  # NLP

    m = Model(Ipopt.Optimizer)

    build_model!(m,p)

    ij_edges = [string(i*"-"*j) for j in p.busses for i in i_to_j(j, p)];

    @objective(m, Min, 
        sum( m[:lij][i_j,t] for t in 1:p.Ntimesteps, i_j in  ij_edges)
    )

    optimize!(m)
    
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

    vs = get_bus_values(:vsqrd, m, p)
    
    for b in keys(vs)
        @test abs(vs[b][1] - dss_voltages[b][1]) < 0.01
    end
    #=
    All BFM vs are slightly higher, which could be explained by not modeling shunts
    =#

end

end  # all tests
