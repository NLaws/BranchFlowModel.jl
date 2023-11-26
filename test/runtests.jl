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
using Graphs, MetaGraphs  # TODO export what is needed from MetaGraphs in BranchFlowModel

# # hack for local testing
# using Pkg
# Pkg.activate("..")
# using BranchFlowModel
# Pkg.activate(".")

Random.seed!(42)


function dss_voltages_pu()
    d = Dict()
    for b in OpenDSSDirect.Circuit.AllBusNames() 
        OpenDSSDirect.Circuit.SetActiveBus(b)
        d[b] = OpenDSSDirect.Bus.puVmagAngle()[1:2:end]
    end
    return d
end


function build_min_loss_model(p)
    m = Model(Ipopt.Optimizer)
    build_model!(m,p)
    @objective(m, Min, 
        sum( m[:lij][i_j,t] for t in 1:p.Ntimesteps, i_j in  p.edge_keys)
    )
    set_optimizer_attribute(m, "print_level", 0)
    return m
end


@testset "BranchFlowModel.jl" begin


@testset "CommonOPF.Network" begin

    fp = joinpath("data", "ieee13", "ieee13_single_phase.yaml")
    net = Network(fp)
    net.v0 = 1.0435118162902168
    m = Model(Ipopt.Optimizer)
    build_model!(m, net; relaxed=false)
    @objective(m, Min, 
        sum( m[:lij][edge,t] for t in 1:net.Ntimesteps, edge in edges(net))
    )
    optimize!(m)

    vs_net = get_variable_values(:vsqrd, m, net)
    lij_net = get_variable_values(:lij, m, net)

    p = Inputs(
        joinpath("data", "ieee13_makePosSeq", "Master.dss"), 
        "650";
        Sbase=net.Sbase, 
        Vbase=net.Vbase, 
        v0 = 1.0,
        v_uplim = 1.05,
        v_lolim = 0.95,
        Ntimesteps = 1,
        relaxed = false  # NLP
    );
    p.regulators[("650", "rg60")][:turn_ratio] = 1.0435118162902168  # dss_voltages["rg60"][1]
    m = build_min_loss_model(p)
    optimize!(m)
    
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

    vs = get_variable_values(:vsqrd, m, p)
    lij = get_variable_values(:lij, m, p)
    
    for b in keys(vs)
        if b in keys(vs_net)
            if b == "650" continue end  
            # ignore this value b/c Inputs model has an extra line 650-rg60 from a transformer
            @test abs(vs[b][1] - vs_net[b][1]) < 0.0001
        end
    end

    for edge in keys(lij_net)
        if edge in keys(lij)
            @test abs(lij[string(edge[1]*"-"*edge[2])][1] - lij_net[edge][1]) < 0.0001
        end
    end


    # TODO use https://jso.dev/NLPModelsJuMP.jl/dev/tutorial/ to newton-solve the NLP with no
    # relaxation and compare

end


@testset "Papavasiliou 2018 with shunts" begin
    #=
    Copied network from paper "Analysis of DLMPs"
    and testing some DLMP values here as well as the addition of shunt susceptance values
    =#
    
    T = 1
    # edges_as_drawn = [
    #     ("0", "1"), ("1", "2"), ("2", "3"), ("3", "4"), ("4", "5"), ("5", "6"),
    #     ("8", "7"), ("3", "8"), 
    #     ("8", "9"), ("9", "10"), ("10", "11"), 
    #     ("0", "12"), ("12", "13"), ("13", "14"), 
    # ]
    # edges sorted according to ending bus
    # and 7 switched with 8, gives closer voltages
    edges = [
        ("0", "1"), ("1", "2"), ("2", "3"), ("3", "4"), ("4", "5"), ("5", "6"),
        ("3", "7"), ("7", "8"), 
        ("7", "9"), ("9", "10"), ("10", "11"), 
        ("0", "12"), ("12", "13"), ("13", "14"), 
    ]
    linecodes = ["l"*string(i) for i = 1:14]
    linelengths = repeat([1.0], length(edges))
    phases = repeat([[1]], length(edges))
    substation_bus = "0"

    Pload = Dict(
        "1" => [0.7936], 
        "2" => [0.0], 
        "3" => [0.0201], 
        "4" => [0.0173], 
        "5" => [0.0291], 
        "6" => [0.0219], 
        "7" => [-0.1969], # uncontrolled generator
        "8" => [0.0235], 
        "9" => [0.0229], 
        "10" => [0.0217], 
        "11" => [0.0132], 
        "12" => [0.6219], 
        "13" => [0.0014], 
        "14" => [0.0224], 
    )

    Qload = Dict(
        "1" => [0.01855], 
        "2" => [0.0], 
        "3" => [0.0084], 
        "4" => [0.0043], 
        "5" => [0.0073], 
        "6" => [0.0055], 
        "7" => [0.0019], 
        "8" => [0.0059], 
        "9" => [0.0142], 
        "10" => [0.0065], 
        "11" => [0.0033], 
        "12" => [0.1291], 
        "13" => [0.0008], 
        "14" => [0.0083], 
    )

    Zdict = Dict(
        "l1" => Dict("rmatrix"=> [0.001], "xmatrix"=> [0.12], "nphases"=> 1),
        "l2" => Dict("rmatrix"=> [0.0883], "xmatrix"=> [0.1262], "nphases"=> 1),
        "l3" => Dict("rmatrix"=> [0.1384], "xmatrix"=> [0.1978], "nphases"=> 1),
        "l4" => Dict("rmatrix"=> [0.0191], "xmatrix"=> [0.0273], "nphases"=> 1),
        "l5" => Dict("rmatrix"=> [0.0175], "xmatrix"=> [0.0251], "nphases"=> 1),
        "l6" => Dict("rmatrix"=> [0.0482], "xmatrix"=> [0.0689], "nphases"=> 1),
        "l7" => Dict("rmatrix"=> [0.0523], "xmatrix"=> [0.0747], "nphases"=> 1),
        "l8" => Dict("rmatrix"=> [0.0407], "xmatrix"=> [0.0582], "nphases"=> 1),
        "l9" => Dict("rmatrix"=> [0.01], "xmatrix"=> [0.0143], "nphases"=> 1),
        "l10" => Dict("rmatrix"=> [0.0241], "xmatrix"=> [0.0345], "nphases"=> 1),
        "l11" => Dict("rmatrix"=> [0.0103], "xmatrix"=> [0.0148], "nphases"=> 1),
        "l12" => Dict("rmatrix"=> [0.001], "xmatrix"=> [0.12], "nphases"=> 1),
        "l13" => Dict("rmatrix"=> [0.1559], "xmatrix"=> [0.1119], "nphases"=> 1),
        "l14" => Dict("rmatrix"=> [0.0953], "xmatrix"=> [0.0684], "nphases"=> 1),
    )
    v0 = 1.0
    shunts = Dict(
        "0" => 0.0,
        "1" => 1.1,
        "2" => 2.8,
        "3" => 2.4,
        "4" => 0.4,
        "5" => 0.8,
        "6" => 0.6,
        "7" => 0.6,
        "8" => 1.2,
        "9" => 0.4,
        "10" => 0.4,
        "11" => 0.1,
        "12" => 0.1,
        "13" => 0.2,
        "14" => 0.1,
    )
    shunts = Dict(k => v*1e-3 for (k,v) in shunts)

    p = Inputs(
        edges, 
        linecodes, 
        linelengths, 
        phases,
        substation_bus;
        Pload=Pload, 
        Qload=Qload, 
        Sbase=1,
        Vbase=1,
        Zdict=Zdict, 
        v0=v0, 
        v_lolim=0.9, 
        v_uplim=1.1,
        Ntimesteps=T, 
        P_up_bound=1e4,
        Q_up_bound=1e4,
        P_lo_bound=-1e4,
        Q_lo_bound=-1e4,
        Isquared_up_bounds=Dict{String, Float64}(),
        relaxed=true,  # TODO does the unrelaxed model match unrelaxed
        shunt_susceptance=shunts,
    );
    # TODO check_soc_inequalities
    # TODO line limits
    m = Model(ECOS.Optimizer) 
    set_optimizer_attribute(m, "maxit", 10000)
    set_optimizer_attribute(m, "verbose", 0)
    set_optimizer_attribute(m, "reltol", 1e-7)
    build_model!(m,p)

    # put in the generator at bus 11
    b = "11"
    @variable(m, 0.4 >= pgen11 >= 0)
    @variable(m, 0.4 >= qgen11 >= 0)
    JuMP.delete.(m, m[:injectioncons][b]["p"])
    m[:injectioncons][b]["p"] = @constraint(m, 
        m[:Pj][b, 1] == pgen11 - p.Pload[b][1]
    )
    JuMP.delete.(m, m[:injectioncons][b]["q"])
    m[:injectioncons][b]["q"] = @constraint(m, 
        m[:Qj][b, 1] == qgen11 - p.Qload[b][1]
    )

    @objective(m, Min,
         m[:Pj][p.substation_bus, 1] * 50 + pgen11 * 10
    )
    optimize!(m)
    r = Results(m,p)

    r.real_power_injections["0"]
    #   1.064072, looking for 1.063

    #=
        Increase load on 11 s.t. generator can't export
        and test if all prices > 50
    =#
    m = Model(ECOS.Optimizer)
    set_optimizer_attribute(m, "maxit", 100000)
    set_optimizer_attribute(m, "verbose", 0)
    set_optimizer_attribute(m, "reltol", 1e-6)

    p.Pload["11"][1]  = 0.4 # this is not enough to get prices > 50 b/c 7 is injecting
    p.Pload["7"][1] = 0  # so set load at 7 to zero

    build_model!(m,p)

    # put in the generator at bus 11
    b = "11"
    @variable(m, 0.4 >= pgen11 >= 0)
    @variable(m, 0.4 >= qgen11 >= 0)
    JuMP.delete.(m, m[:injectioncons][b]["p"])
    m[:injectioncons][b]["p"] = @constraint(m, 
        m[:Pj][b, 1] == pgen11 - p.Pload[b][1]
    )
    JuMP.delete.(m, m[:injectioncons][b]["q"])
    m[:injectioncons][b]["q"] = @constraint(m, 
        m[:Qj][b, 1] == qgen11 - p.Qload[b][1]
    )

    @objective(m, Min,
         m[:Pj][p.substation_bus, 1] * 50 + pgen11 * 10
    )
    optimize!(m)
    r = Results(m,p)
    @test all((price[1] >= 49.999 for price in values(r.shadow_prices)))

    #=
     now with voltage limit preventing full gen the gen should set the price
    =#
    p.v_uplim = 1.01 # pgen11 not curtailed at 1.03 ?!
    p.Pload["11"][1] = 0.0132  # back to original value
    p.relaxed = false
    m = Model(Ipopt.Optimizer)
    set_optimizer_attribute(m, "print_level", 0)

    build_model!(m,p)

    # put in the generator at bus 11
    b = "11"
    @variable(m, 0.4 >= pgen11 >= 0)
    @variable(m, 0.4 >= qgen11 >= 0)
    JuMP.delete.(m, m[:injectioncons][b]["p"])
    m[:injectioncons][b]["p"] = @constraint(m, 
        m[:Pj][b, 1] == pgen11 - p.Pload[b][1]
    )
    JuMP.delete.(m, m[:injectioncons][b]["q"])
    m[:injectioncons][b]["q"] = @constraint(m, 
        m[:Qj][b, 1] == qgen11 - p.Qload[b][1]
    )

    @objective(m, Min,
         m[:Pj][p.substation_bus, 1] * 50 + pgen11 * 10
    )
    optimize!(m)
    r = Results(m,p)
    @test r.shadow_prices["11"][1] ≈ 10.0


    #=
    is the DLMP lower with DER not net injecting?
    i.e. is the total DSO cost lower when paying DLMP > LMP?
    =#

    p.v_uplim = 1.05 # pgen11 not curtailed at 1.03 ?!
    p.Pload["11"][1] = 0.1 
    p.relaxed = false
    m = Model(Ipopt.Optimizer)
    set_optimizer_attribute(m, "print_level", 0)

    build_model!(m,p)
    @objective(m, Min,
         m[:Pj][p.substation_bus, 1] * 50
    )
    optimize!(m)
    prices_no_der = Dict(
        j => JuMP.dual.(m[:loadbalcons][j]["p"][1])
        for j in p.busses
    )  # TODO put this as option in Results?
    cost_no_der = objective_value(m)


    m = Model(Ipopt.Optimizer)
    set_optimizer_attribute(m, "print_level", 0)
    build_model!(m,p)
    # put in the generator at bus 11
    b = "11"
    @variable(m, 0.05 >= pgen11 >= 0)
    @variable(m, 0.05 >= qgen11 >= 0)
    JuMP.delete.(m, m[:injectioncons][b]["p"])
    m[:injectioncons][b]["p"] = @constraint(m, 
        m[:Pj][b, 1] == pgen11 - p.Pload[b][1]
    )
    JuMP.delete.(m, m[:injectioncons][b]["q"])
    m[:injectioncons][b]["q"] = @constraint(m, 
        m[:Qj][b, 1] == qgen11 - p.Qload[b][1]
    )

    @objective(m, Min,
         m[:Pj][p.substation_bus, 1] * 50 + pgen11 * 10
    )
    optimize!(m)
    r = Results(m,p)
    cost_with_der = value(m[:Pj][p.substation_bus, 1]) * 50 + r.shadow_prices["11"][1] * value(pgen11)
    @test cost_with_der < cost_no_der
    @test r.shadow_prices["11"][1] < prices_no_der["11"][1]

end


@testset "Results" begin
    Sbase = 1e6
    Vbase = 12.47e3
    p = Inputs(
        joinpath("data", "singlephase38lines", "master.dss"), 
        "0";
        Sbase=Sbase, 
        Vbase=Vbase, 
        v0 = 1.00,
        v_uplim = 1.05,
        v_lolim = 0.95,
        relaxed = false,
    );
    m = build_min_loss_model(p)
    optimize!(m)
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]
    r = Results(m, p)
end


@testset "merge parallel single phase lines" begin
    #=       3
           c -- e                   
         2/      1\                 6.5  
    a -- b         g    ->   a -- b -- g
       2.5\      2/                     
           d -- f            
            3.5
    Merge parallel lines sets that do not have loads
    =#
    
    edges = [("a", "b"), ("b", "c"), ("b", "d"), ("c", "e"), ("d", "f"), ("e", "g"), ("f", "g")]
    linecodes = ["l1", "l2", "l3", "l2", "l3", "l2", "l3"]
    linelengths = [1.0, 2.0, 2.5, 3.0, 3.5, 1.0, 1.0]
    phases = [[1,2], [1], [2], [1], [2], [1], [1,2]]
    substation_bus = "a"
    Pload = Dict()
    Qload = Dict()
    Zdict = Dict(
        "l1" => Dict("rmatrix"=> [[1.0, 0.5], [0.5, 1.0]], "xmatrix"=> [[1.0, 0.5], [0.5, 1.0]], "nphases"=> 2),
        "l2" => Dict("rmatrix"=> [2.0], "xmatrix"=> [2.0], "nphases"=> 1),  # total R = 2 * 6
        "l3" => Dict("rmatrix"=> [3.0], "xmatrix"=> [3.0], "nphases"=> 1),  # total R = 3 * 7
    )
    v0 = 1.0
    # NOTE intermediate steps are tested in CommonOPF

    p = Inputs(
        edges, 
        linecodes, 
        linelengths, 
        phases,
        substation_bus;
        Pload=Pload, 
        Qload=Qload, 
        Sbase=1, 
        Vbase=1, 
        Zdict=Zdict, 
        v0=v0, 
        Isquared_up_bounds=Dict{String, Float64}()
    )

    combine_parallel_lines!(p)
    @test p.busses == ["a", "b", "g"]
    @test p.edge_keys == ["a-b", "b-g"]
    @test get_ijlinelength("b", "g", p) == 6.5  # avg of 6 and 7
    @test rij("b", "g", p)[1,1] == 18
    @test rij("b", "g", p)[2,2] == 21
    @test zij("b", "g", p)[1,1] == 18+18im
    @test zij("b", "g", p)[2,2] == 21+21im

end
    


@testset "ieee13 positive sequence" begin
    # make the dss solution to compare
    dss("Redirect data/ieee13/IEEE13Nodeckt.dss")
    @test(OpenDSSDirect.Solution.Converged() == true)
    dss_voltages = dss_voltages_pu()

    p3phase = Inputs(joinpath("data", "ieee13", "IEEE13Nodeckt.dss"), "rg60");

    extract_phase = 3
    p = Inputs(
        joinpath("data", "ieee13", "IEEE13Nodeckt.dss"), 
        "rg60";
        Sbase = 1_000_000, 
        Vbase = 4160/sqrt(3), 
        v0 = dss_voltages["rg60"][extract_phase],
        v_uplim = dss_voltages["rg60"][extract_phase],
        v_lolim = minimum(values(dss_voltages))[1],
        Ntimesteps = 1,
        extract_phase = extract_phase,
        relaxed = false
    );

    @test check_connected_graph(p) == true

    m = build_min_loss_model(p)
    optimize!(m)
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

    # vs = get_variable_values(:vsqrd, m, p)

    # for (bus, phs) in p3phase.phases_into_bus
    #     if bus in keys(vs)
    #         dss_index = indexin(extract_phase, phs)[1]
    #         if !isnothing(dss_index)
    #             println(rpad(bus,5), 
    #                 rpad(round(vs[bus][1] -  dss_voltages[bus][dss_index], digits=4), 8)
    #             )
    #         end
    #     end
    # end

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
    @objective(m, Min, sum(m[:Pj]["0", t] for t=1:T))
    optimize!(m)

    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]

end


@testset "ieee13 unbalanced MultiPhase" begin

    # make the dss solution to compare
    dss("Redirect data/ieee13/IEEE13Nodeckt.dss")
    @test(OpenDSSDirect.Solution.Converged() == true)

    dss_voltages = dss_voltages_pu()

    vbase = 4160/sqrt(3)
    p = Inputs(
        joinpath("data", "ieee13", "IEEE13Nodeckt.dss"), 
        "rg60";
        Sbase=1_000_000, 
        Vbase=vbase, 
        v0 = dss_voltages["rg60"],  # openDSS rg60 values
        v_uplim = 1.06,
        v_lolim = 0.94,
        Ntimesteps = 1
    );
    # p.Isquared_up_bounds = Dict(
    #     lc => 100 for lc in Set(p.linecodes)
    # )
    p.P_lo_bound = -10
    p.Q_lo_bound = -10
    p.P_up_bound = 10
    p.Q_up_bound = 10

    m = Model(CSDP.Optimizer)
    set_attribute(m, "printlevel", 0)

    # m = Model(COSMO.Optimizer)

    # m = Model(SCS.Optimizer)

    build_model!(m,p)

    @objective(m, Min, 
        sum( sum(real.(diag(m[:l][t][i_j]))) for t in 1:p.Ntimesteps, i_j in  p.edge_keys)
    )

    optimize!(m)
    
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]

    @test_nowarn(check_rank_one(m,p))

    vs = Dict(k => real.(diag(v)) for (k,v) in voltage_values_by_time_bus(m,p)[1])

    # I = get_variable_values(:l, m, p)  # TODO method for multiphase results

    # for b in keys(vs)
    #     for (i,phsv) in enumerate(filter(v -> v != 0, vs[b]))
    #         # @assert abs(phsv - dss_voltages[b][i]) < 1e-2 "bus $b phase $i failed"
    #         println("$b - $i  $(phsv - dss_voltages[b][i])")
    #     end
    # end
    

    # TODO why are BFM voltages not matching openDSS well?

end


@testset "ieee13 balanced SinglePhase" begin

    # make the dss solution to compare
    dss("clear")
    dss("Redirect data/ieee13_makePosSeq/Master.dss")
    dss("Solve")
    @test(OpenDSSDirect.Solution.Converged() == true)

    dss_voltages = dss_voltages_pu()

    vbase = 4160/sqrt(3)
    p = Inputs(
        joinpath("data", "ieee13_makePosSeq", "Master.dss"), 
        "650";
        Sbase=1_000_000, 
        Vbase=vbase, 
        v0 = 1.0,
        v_uplim = 1.05,
        v_lolim = 0.95,
        Ntimesteps = 1,
        relaxed = false  # NLP
    );
    @test check_connected_graph(p) == true
    g = make_graph(p.busses, p.edges; directed=true)
    @test g.graph.ne == Graphs.nv(g) - 1  # a radial network has n_edges = n_vertices - 1
    @test_warn "The per unit impedance values should be much less than one" check_unique_solution_conditions(p)
    p.regulators[("650", "rg60")][:turn_ratio] = dss_voltages["rg60"][1]

    m = build_min_loss_model(p)
    optimize!(m)
    
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

    vs = get_variable_values(:vsqrd, m, p)
    I = get_variable_values(:lij, m, p)
    
    for b in keys(vs)
        @test abs(vs[b][1] - dss_voltages[b][1]) < 0.01
    end

    bs, depths = busses_from_deepest_to_source(g, "650")
    @test depths[end] == 0 && bs[end] == "650"
    @test depths[end-1] == 1 && bs[end-1] == "rg60"
    @test "675" in bs[1:3] && depths[1] == 6
    @test "652" in bs[1:3] && depths[2] == 6
    @test "611" in bs[1:3] && depths[3] == 6

    splitting_bs, subgraph_bs = splitting_busses(p, "650"; max_busses=7)  # 671 and rg60
    @test length(splitting_bs) == 2
    @test splitting_bs[1] == "671" && splitting_bs[2] == "rg60"
    @test length(subgraph_bs) == 2
    @test length(subgraph_bs[1]) == 7
    @test length(subgraph_bs[2]) == 7
    @test isempty(intersect(subgraph_bs[1], subgraph_bs[2]))

    # the subgraph_bs do not include over laps if add_connections=false
    mg = split_at_busses(p, splitting_bs, subgraph_bs; add_connections=false)
    @test !("671" in mg[3,:p].busses)

    mg = split_at_busses(p, splitting_bs, subgraph_bs; add_connections=true)
    @test length(mg[1,:p].busses) == 2
    @test length(mg[2,:p].busses) == 7
    @test length(mg[3,:p].busses) == 8  # get one more bus then max_busses from adding connection @ 671

    # to get overlaps use connect_subgraphs_at_busses (which is an option in split_at_busses)
    new_subgraphs = connect_subgraphs_at_busses(p, splitting_bs, subgraph_bs)
    mg = split_at_busses(p, splitting_bs, new_subgraphs; add_connections=false)
    @test length(mg[1,:p].busses) == 2
    @test length(mg[2,:p].busses) == 7
    @test length(mg[3,:p].busses) == 8  # get one more bus then max_busses from adding connection @ 671

    @test length(new_subgraphs) == 2
    @test intersect(new_subgraphs[1], new_subgraphs[2])[1] == "671"
    @test intersect(mg[2,:p].busses, mg[3,:p].busses)[1] == "671"
    @test intersect(mg[1,:p].busses, mg[3,:p].busses)[1] == "rg60"

    # test the decomposed solution against openDSS
    builder = Dict(
        v => build_min_loss_model for v in vertices(mg)
    )
    solve_metagraph!(mg, builder, [1e-3, 1e-4, 1e-4]; verbose=false)

    vs = metagraph_voltages(mg)
    for b in keys(vs)
        @test abs(vs[b][1] - dss_voltages[b][1]) < 0.01
    end

end


@testset "SinglePhase network reduction" begin

    function make_solve_min_loss_model(p)
        m = Model(Ipopt.Optimizer)
        build_model!(m,p)
        @objective(m, Min, 
            sum( m[:lij][i_j,t] for t in 1:p.Ntimesteps, i_j in  p.edge_keys)
        )
        set_optimizer_attribute(m, "print_level", 0)
        optimize!(m)
        return m
    end

    # 1 validate BFM against OpenDSS
    Sbase = 1e6
    Vbase = 12.47e3
    p = Inputs(
        joinpath("data", "singlephase38lines", "master.dss"), 
        "0";
        Sbase=Sbase, 
        Vbase=Vbase, 
        v0 = 1.00,
        v_uplim = 1.05,
        v_lolim = 0.95,
        relaxed = false,
    );
    m = make_solve_min_loss_model(p)
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]
    vs = get_variable_values(:vsqrd, m, p)
    # make the dss solution to compare
    dss("clear")
    dss("Redirect data/singlephase38lines/master.dss")
    dss("Solve")
    @test(OpenDSSDirect.Solution.Converged() == true)
    dss_voltages = dss_voltages_pu()
    for b in keys(vs)
        @test abs(vs[b][1] - dss_voltages[b][1]) < 0.001
    end
    nvar_original = JuMP.num_variables(m)

    # 2 validate BFM results stay the same after reduction
    nbusses_before = length(p.busses)
    reduce_tree!(p)
    removed_busses = ["1", "4", "8", "26", "29", "14", "16", "20", "21", "23", "24"]
    for b in removed_busses
        @test !(b in p.busses)
    end
    @test nbusses_before - length(removed_busses) == length(p.busses)
    m = make_solve_min_loss_model(p)
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]
    vs_reduced = get_variable_values(:vsqrd, m, p)
    for b in keys(vs_reduced)
        @test abs(dss_voltages[b][1] - vs_reduced[b][1]) < 0.001
    end
    nvar_reduced = JuMP.num_variables(m)
    @test nvar_original > nvar_reduced  # 225 > 159

    # 3 split and solve the reduced model, compare v
    g = BranchFlowModel.make_graph(p.busses, p.edges; directed=true)
    p_above, p_below = BranchFlowModel.split_inputs(p, "12");
    @test intersect(p_above.busses, p_below.busses) == ["12"]
    @test length(p.busses) == length(p_below.busses) + length(p_above.busses) - 1
    @test isempty(intersect(p_above.edges, p_below.edges))
    @test length(p.edges) == length(p_below.edges) + length(p_above.edges)
    # splitting at 12 should put 12-25 in p_below (except for the removed busses)
    for b in setdiff!(string.(12:25), removed_busses)
        @test b in p_below.busses
    end
    # solve above first with sum of p_below loads, then set p_below.v0, solve p_below, set p_above.P/Qload to p_below.substation_bus values
    # later solve in parallel
    # p_above.Pload and Qload exist already for bus "12"; we need to add p_below's loads to it
    init_inputs!([p_below, p_above])  # set p_above loads at p_below.substation_bus
    # NOTE only one time step so all load vectors have one value
    @test p_above.Pload["12"][1] == sum( v[1] for v in values(p_below.Pload) )
    @test p_above.Qload["12"][1] == sum( v[1] for v in values(p_below.Qload) )

    m_above = make_solve_min_loss_model(p_above)
    init_vs = Dict(
        p_below.substation_bus => sqrt(value(m_above[:vsqrd][p_below.substation_bus,1]))
    )
    init_inputs!([p_above, p_below]; init_vs=init_vs)
    @test p_below.v0 == sqrt(value(m_above[:vsqrd][p_below.substation_bus,1]))

    m_below = make_solve_min_loss_model(p_below)

    # at this point the v's at bus 12 agree by design
    # but the loads do not b/c we did not account for losses in first iteration
    # the plus sign is not a minus because the P/Q values should be equal and _opposite_
    pdiff = value(m_below[:Pj]["12",1]) + value(m_above[:Pj]["12",1])  # 0.001606
    qdiff = value(m_below[:Qj]["12",1]) + value(m_above[:Qj]["12",1])  # 0.000526
    vdiff = 1.0

    tol = 1e-6
    while pdiff > tol || qdiff > tol || vdiff > tol
        p_above.Pload["12"][1] = value(m_below[:Pj]["12",1]) * p_above.Sbase
        p_above.Qload["12"][1] = value(m_below[:Qj]["12",1]) * p_above.Sbase
        m_above = make_solve_min_loss_model(p_above)
        v_above = sqrt(value(m_above[:vsqrd][p_below.substation_bus,1]))
        vdiff = abs(p_below.v0 - v_above)
        p_below.v0 = v_above
        m_below = make_solve_min_loss_model(p_below)
        pdiff = value(m_below[:Pj]["12",1]) + value(m_above[:Pj]["12",1])  # 3.758466e-8
        qdiff = value(m_below[:Qj]["12",1]) + value(m_above[:Qj]["12",1])  # 1.231675e-8
    end

    vs_decomposed = get_variable_values(:vsqrd, m_above, p_above)
    merge!(vs_decomposed, get_variable_values(:vsqrd, m_below, p_below))
    for b in keys(vs_decomposed)
        @test abs(dss_voltages[b][1] - vs_decomposed[b][1]) < 0.001
    end

    # split model into 3 models and solve
    mg = split_at_busses(p, ["7", "13"])
    @test mg[1, :p].substation_bus == "0"
    @test mg[2, :p].substation_bus == "7"
    @test mg[3, :p].substation_bus == "13"
    # init_inputs!(mg)  # done in split_at_busses
    for v in get_prop(mg, :load_sum_order)
        set_prop!(mg, v, :m, make_solve_min_loss_model(mg[v, :p]))
    end
    set_indexing_prop!(mg, :m)
    # every model solved once using load approximations
    # now set better load approximations and connect voltages
    set_inputs!(mg)
    @test mg[2,:p].v0[1] ≈ sqrt(value(mg[1,:m][:vsqrd][mg[2,:p].substation_bus,1]))
    @test mg[3,:p].v0[1] ≈ sqrt(value(mg[2,:m][:vsqrd][mg[3,:p].substation_bus,1]))
    # run another solve
    for v in get_prop(mg, :load_sum_order)
        set_prop!(mg, v, :m, make_solve_min_loss_model(mg[v, :p]))
    end
    pdiffs1, qdiffs1, vdiffs1 = get_diffs(mg)
    # another round to compare:
    set_inputs!(mg)
    for v in get_prop(mg, :load_sum_order)
        set_prop!(mg, v, :m, make_solve_min_loss_model(mg[v, :p]))
    end
    pdiffs2, qdiffs2, vdiffs2 = get_diffs(mg)
    @test sum(pdiffs2) < sum(pdiffs1) 
    @test sum(qdiffs2) < sum(qdiffs1) 
    @test sum(vdiffs2) < sum(vdiffs1) 

    vs_decomposed = get_variable_values(:vsqrd, mg[1, :m], mg[1, :p])
    merge!(vs_decomposed, get_variable_values(:vsqrd, mg[2, :m], mg[2, :p]))
    merge!(vs_decomposed, get_variable_values(:vsqrd, mg[3, :m], mg[3, :p]))
    for b in keys(vs_decomposed)
        @test abs(dss_voltages[b][1] - vs_decomposed[b][1]) < 0.001
    end
end


@testset "Single phase regulators at network splits" begin
    Sbase = 1e6
    Vbase = 12.47e3
    # use OpenDSSDirect to replace a line 20-21 with a transformer/regulator and test solution convergence
    p = Inputs(
        joinpath("data", "singlephase38lines", "master_extra_trfx.dss"), 
        "0";
        Sbase=Sbase, 
        Vbase=Vbase, 
        v0 = 1.00,
        v_uplim = 1.05,
        v_lolim = 0.95,
        relaxed = false,
    );
    @test ("20","21") in keys(p.regulators)
    @test reg_busses(p) == ["21"]

    dss("clear")
    dss("Redirect data/singlephase38lines/master_extra_trfx.dss")
    dss("Solve")
    @test(OpenDSSDirect.Solution.Converged() == true)
    dss_voltages = dss_voltages_pu()

    p.regulators[("20","21")][:turn_ratio] = dss_voltages["21"][1] / dss_voltages["20"][1]
    @test turn_ratio(p, "21") == dss_voltages["21"][1] / dss_voltages["20"][1]

    # split at two busses including the regulator and compare against openDSS
    mg = split_at_busses(p, ["14","21"])
    builder = Dict(
        v => build_min_loss_model for v in vertices(mg)
    )
    solve_metagraph!(mg, builder, [1e-5, 1e-5, 1e-5]; verbose=false)
    vs = metagraph_voltages(mg)
    for b in keys(vs)
        @test abs(vs[b][1] - dss_voltages[b][1]) < 0.001
    end

    # split only at regulator, add vreg, and test the vdiff is zero
    # (because using vreg makes the regulator voltage equal to the setting)
    p.regulators[("20","21")][:vreg] = 1.02
    mg = split_at_busses(p, ["21"])

    builder = Dict(
        v => build_min_loss_model for v in vertices(mg)
    )
    solve_metagraph!(mg, builder, [1e-3, 1e-3, 1e-3]; verbose=false)
    pdiffs, qdiffs, vdiffs = get_diffs(mg)
    @test maximum(vdiffs) ≈ 0  # because that is how it is defined across a regulator


end


@testset "split_at_busses" begin
    #=     c -- e                    c | c -- e
          /                         /
    a -- b           ->   a -- b | b
          \                         \
           d -- f                    d | d -- f
    =#

    edges = [("a", "b"), ("b", "c"), ("b", "d"), ("c", "e"), ("d", "f")]
    linecodes = repeat(["l1"], length(edges))
    linelengths = repeat([1.0], length(edges))
    phases = repeat([[1]], length(edges))
    substation_bus = "a"
    Pload = Dict("c" => [1.0], "d" => [1.0], "e" => [1.0], "f" => [1.0])
    Qload = Dict("c" => [0.1], "d" => [0.1], "e" => [0.1], "f" => [0.1])
    Zdict = Dict("l1" => Dict("rmatrix"=> [1.0], "xmatrix"=> [1.0], "nphases"=> 1))
    v0 = 1.0

    p = Inputs(
        edges, 
        linecodes, 
        linelengths, 
        phases,
        substation_bus;
        Pload=Pload, 
        Qload=Qload, 
        Sbase=1, 
        Vbase=1, 
        Zdict=Zdict, 
        v0=v0, 
        v_lolim=0.95, 
        v_uplim=1.05,
        Ntimesteps=1, 
        P_up_bound=1e4,
        Q_up_bound=1e4,
        P_lo_bound=-1e4,
        Q_lo_bound=-1e4,
        Isquared_up_bounds=Dict{String, Float64}(),
        relaxed=true
    )

    mg = split_at_busses(p, ["c", "d", "b"])
    @test MetaGraphs.outneighbors(mg, 1) == [4]
    @test MetaGraphs.outneighbors(mg, 4) == [2, 3]
    @test MetaGraphs.outneighbors(mg, 2) == MetaGraphs.outneighbors(mg, 3) == Int[]

    p_above, p_below = split_inputs(p, "b", ["b", "c", "e"])
    @test length(p_below.busses) == 3
    @test "b" in p_below.busses && "c" in p_below.busses && "e" in p_below.busses
    @test length(p_above.busses) == 4
    @test "b" in p_above.busses && "a" in p_above.busses 
    @test "d" in p_above.busses && "f" in p_above.busses
    @test length(p_below.edges) == 2
    @test ("b", "c") in p_below.edges && ("c", "e") in p_below.edges
    @test length(p_above.edges) == 3
    @test ("a", "b") in p_above.edges
    @test ("b", "d") in p_above.edges
    @test ("d", "f") in p_above.edges



    # TODO test CommonOPF.trim_tree!
    delete!(p.Pload, "e")
    delete!(p.Qload, "e")
    BranchFlowModel.CommonOPF.trim_tree!(p)
end

end  # all tests
