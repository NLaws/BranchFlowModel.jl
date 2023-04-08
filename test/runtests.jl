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
using Graphs, MetaGraphs  # maybe export what is needed from MetaGraphs in BranchFlowModel ?

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
    g = make_graph(p.busses, p.edges)
    @test g.graph.ne == length(g.vprops) - 1  # a radial network has n_edges = n_vertices - 1
    @test_warn "The per unit impedance values should be much less than one" check_unique_solution_conditions(p)
    p.regulators[("650", "rg60")][:turn_ratio] = dss_voltages["rg60"][1]

    m = build_min_loss_model(p)
    optimize!(m)
    
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

    vs = get_bus_values(:vsqrd, m, p)
    
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
    vs = get_bus_values(:vsqrd, m, p)
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
    vs_reduced = get_bus_values(:vsqrd, m, p)
    for b in keys(vs_reduced)
        @test abs(dss_voltages[b][1] - vs_reduced[b][1]) < 0.001
    end
    nvar_reduced = JuMP.num_variables(m)
    @test nvar_original > nvar_reduced  # 225 > 159

    # 3 split and solve the reduced model, compare v
    g = BranchFlowModel.make_graph(p.busses, p.edges)
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

    vs_decomposed = get_bus_values(:vsqrd, m_above, p_above)
    merge!(vs_decomposed, get_bus_values(:vsqrd, m_below, p_below))
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

    vs_decomposed = get_bus_values(:vsqrd, mg[1, :m], mg[1, :p])
    merge!(vs_decomposed, get_bus_values(:vsqrd, mg[2, :m], mg[2, :p]))
    merge!(vs_decomposed, get_bus_values(:vsqrd, mg[3, :m], mg[3, :p]))
    for b in keys(vs_decomposed)
        @test abs(dss_voltages[b][1] - vs_decomposed[b][1]) < 0.001
    end

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

    dss("clear")
    dss("Redirect data/singlephase38lines/master_extra_trfx.dss")
    dss("Solve")
    @test(OpenDSSDirect.Solution.Converged() == true)
    dss_voltages = dss_voltages_pu()

    p.regulators[("20","21")][:turn_ratio] = dss_voltages["21"][1] / dss_voltages["20"][1]

    mg = split_at_busses(p, ["21"])

    builder = Dict(
        v => build_min_loss_model for v in vertices(mg)
    )
    solve_metagraph!(mg, builder, [1e-3, 1e-3, 1e-3]; verbose=false)
    pdiffs, qdiffs, vdiffs = get_diffs(mg)
    @test maximum(vdiffs) ≈ 0  # because that is how it is defined across a regulator

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
    Zdict = Dict("l1" => Dict("rmatrix"=> [1.0], "zmatrix"=> [1.0], "nphases"=> 1))
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



    # TODO test trim_tree!
    delete!(p.Pload, "e")
    delete!(p.Qload, "e")
    trim_tree!(p)
end

end  # all tests
