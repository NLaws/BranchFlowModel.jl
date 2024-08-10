using BranchFlowModel
using Test
using Random
using ECOS
using JuMP
import OpenDSSDirect as OpenDSS
using SCS
using LinearAlgebra
using COSMO
using CSDP
using Ipopt
import Graphs

# # hack for local testing
# using Pkg
# Pkg.activate("..")
# using BranchFlowModel
# Pkg.activate(".")

CPF = BranchFlowModel.CommonOPF

Random.seed!(42)


function dss_voltages_pu()
    d = Dict()
    for b in OpenDSS.Circuit.AllBusNames() 
        OpenDSS.Circuit.SetActiveBus(b)
        d[b] = OpenDSS.Bus.puVmagAngle()[1:2:end]
    end
    return d
end



function dss_voltages_mag_angle()
    d = Dict()
    for b in OpenDSS.Circuit.AllBusNames() 
        OpenDSS.Circuit.SetActiveBus(b)
        d[b] = OpenDSS.Bus.VMagAngle()#[1:2:end]
    end
    return d
end


function build_min_loss_model(net::CPF.Network{CPF.SinglePhase})
    m = Model(Ipopt.Optimizer)
    BranchFlowModel.build_bfm!(m, net; relaxed=false)
    @objective(m, Min, 
        sum( 
            m[:lij][i_j][t] for t in 1:net.Ntimesteps, i_j in CPF.edges(net)
        )
    )
    set_optimizer_attribute(m, "print_level", 0)
    return m
end


function hermitian_variable_to_vector(m::JuMP.AbstractModel, var::Symbol, t::Int, bus::Union{String, Tuple{String, String}}; tol=1e-5)

    H = value.(m[var][t][bus])

    # Perform eigenvalue decomposition
    eig = eigen(H)

    # # Eigenvalues and eigenvectors
    eigenvalues = eig.values
    eigenvectors = eig.vectors

    # Select a non-zero eigenvalue and corresponding eigenvector
    non_zero_indices = findall(x -> abs(x) > tol, eigenvalues)
    lambda_i = abs(eigenvalues[non_zero_indices[1]])
    u_i = eigenvectors[:, non_zero_indices[1]]

    # Construct the vector v
    v = sqrt(lambda_i) * u_i
    return v
end



function check_opendss_powers(;tol=1e-6)

    elements = OpenDSS.Circuit.AllElementNames()

    for element in elements
        # Set the active element
        if !startswith(element, "Load")
            continue
        end
        OpenDSS.Circuit.SetActiveElement(element)
        
        # Get the bus name where the element is connected
        bus_name = OpenDSS.CktElement.BusNames()[1]  # Get the first bus name
        
        # Get the power in kW and kvar
        power = OpenDSS.CktElement.TotalPowers()  # Returns power in kW and kvar for each phase

        load_name = string(split(element, ".")[2])
        OpenDSS.Loads.Name(load_name)
        
        p_mismatch = real(power)[1] - OpenDSS.Loads.kW()
        q_mismatch = imag(power)[1] - OpenDSS.Loads.kvar()

        if abs(p_mismatch) > tol
            return false
        end
        if abs(q_mismatch) > tol
            return false
        end
    end

    return true
end


@testset "BranchFlowModel.jl" begin


@testset "utilities" begin

    # matrix_phases_to_vec used to build multiphase KVL
    M = [1 2 3; 4 5 6; 7 8 9]

    phases = [2, 3]
    v = BranchFlowModel.matrix_phases_to_vec(M, phases)
    @test v == [5, 6, 8, 9]

    phases = [1]
    v = BranchFlowModel.matrix_phases_to_vec(M, phases)
    @test v == [1]

    phases = [1, 2, 3]
    v = BranchFlowModel.matrix_phases_to_vec(M, phases)
    @test v == [1, 2, 3, 4, 5, 6, 7, 8, 9]
end


include("test_nlp.jl")


@testset "multiphase KVL" begin

    # simple min loss model
    net = CPF.Network(joinpath("data", "two_line_multi_phase.yaml"))
    net.bounds.i_upper_mag = 15 * net.Sbase / net.Vbase

    m = Model(CSDP.Optimizer)
    set_attribute(m, "printlevel", 0)
    build_bfm!(m, net)
    @objective(m, Min, 
        sum( sum(real.(diag(m[:l][t][i_j]))) for t in 1:net.Ntimesteps, i_j in  edges(net))
    )
    
    optimize!(m)
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]

    w = Dict(
        k => JuMP.value.(w)
        for (k,w) in m[:w][1]
    )

    Lijs = Dict(
        k => JuMP.value.(lij)
        for (k, lij) in m[:l][1]
    )

    Sijs = Dict(
        k => JuMP.value.(Sij)
        for (k, Sij) in m[:Sij][1]
    )

    i = "b1"
    j = "b2"
    @test CPF.phases_into_bus(net, j) == [1,2,3]

    Z = CPF.zij_per_unit(i, j, net)
    Sij = Sijs[(i,j)]
    Lij = Lijs[(i,j)]

    RHS = w[i] - Sij * conj(transpose(Z)) - Z * conj(transpose(Sij)) + Z * Lij * conj(transpose(Z))

    for (a,b) in zip(w[j], RHS)
        @test a ≈ b atol=1e-5
    end

    i = "b2"
    j = "b4"
    phases = CPF.phases_into_bus(net, j)
    @test phases == [1,3]
    T = BranchFlowModel.matrix_phases_to_vec

    Z = CPF.zij_per_unit(i, j, net)
    Sij = Sijs[(i,j)]
    Lij = Lijs[(i,j)]

    RHS = w[i] - Sij * conj(transpose(Z)) - Z * conj(transpose(Sij)) + Z * Lij * conj(transpose(Z))

    for (a,b) in zip( T(w[j], phases), T(RHS, phases) )
        @test a ≈ b atol=1e-5
    end

end


@testset "Papavasiliou 2018 with shunts" begin
    #=
    Copied network from paper "Analysis of DLMPs"
    and testing some DLMP values here as well as the addition of shunt susceptance values
    =#
    
    # edges_as_drawn = [
    #     ("0", "1"), ("1", "2"), ("2", "3"), ("3", "4"), ("4", "5"), ("5", "6"),
    #     ("8", "7"), ("3", "8"), 
    #     ("8", "9"), ("9", "10"), ("10", "11"), 
    #     ("0", "12"), ("12", "13"), ("13", "14"), 
    # ]
    # edges sorted according to ending bus
    # and 7 switched with 8, gives closer voltages

    net = Network_Papavasiliou_2018()
    # TODO check_soc_inequalities
    # TODO line limits

    m_net = Model(ECOS.Optimizer) 
    set_optimizer_attribute(m_net, "maxit", 10000)
    set_optimizer_attribute(m_net, "verbose", 0)
    set_optimizer_attribute(m_net, "reltol", 1e-7)
    build_bfm!(m_net, net)

    function add_generator_at_bus_11_net!(m)
        b = "11"
        @variable(m, 0.4 >= pgen11 >= 0)
        @variable(m, 0.4 >= qgen11 >= 0)
        JuMP.delete.(m, m[:loadbalcons][b]["p"])
        m[:loadbalcons][b]["p"] = @constraint(m, 
            sum( m[:Pij][(i,b)][1] for i in i_to_j(b, net) )
            - sum( m[:lij][(i,b)][1] * rij(i,b,net) for i in i_to_j(b, net) ) 
            + pgen11 - net[b][:Load].kws1[1] == 0
        )
        JuMP.delete.(m, m[:loadbalcons][b]["q"])
        m[:loadbalcons][b]["q"] = @constraint(m, 
            sum( m[:Pij][(i,b)][1] for i in i_to_j(b, net) )
            - sum( m[:lij][(i,b)][1] * rij(i,b,net) for i in i_to_j(b, net) ) 
            + qgen11 - net[b][:Load].kvars1[1] == 0
        )
    end
    add_generator_at_bus_11_net!(m_net)

    @objective(m_net, Min,
         m_net[:p0][1] * 50 + m_net[:pgen11] * 10
    )
    optimize!(m_net)

    vs_net = get_variable_values(:vsqrd, m_net)
    lij_net = get_variable_values(:lij, m_net)
    
    #=
        Increase load on 11 s.t. generator can't export
        and test if all prices > 50
    =#
    m = Model(ECOS.Optimizer)
    set_optimizer_attribute(m, "maxit", 100000)
    set_optimizer_attribute(m, "verbose", 0)
    set_optimizer_attribute(m, "reltol", 1e-6)

    # this is not enough to get prices > 50 b/c 7 is injecting
    net["11"][:Load].kws1 = [0.4]
    # so set load at 7 to zero
    net["7"][:Load].kws1 = [0.0]

    build_bfm!(m, net)

    add_generator_at_bus_11_net!(m)

    @objective(m, Min,
         m[:p0][1] * 50 + m[:pgen11] * 10
    )
    optimize!(m)
    shadow_prices = Dict(
        j => JuMP.dual.(m[:loadbalcons][j]["p"])
        for j in busses(net)
    )

    @test all((price[1] >= 49.999 for price in values(shadow_prices)))

    #=
     with voltage limit preventing full gen the gen should set the price
     TODO translate Papavasiliou power limits to line limits and test
     NOTE have not got the same results as Papavasiliou yet because the network definition is unclear
    =#

end


@testset "ieee13 balanced SinglePhase" begin

    # make the dss solution to compare
    dssfilepath = "data/ieee13_makePosSeq/Master.dss"
    work_dir = pwd()
    OpenDSS.dss("""
        clear
        compile $dssfilepath
        solve
    """)
    cd(work_dir)
    @test(OpenDSS.Solution.Converged() == true)

    dss_voltages = dss_voltages_pu()

    net = CPF.dss_to_Network(dssfilepath)
    @test net.substation_bus == "sourcebus"
    net.Vbase = 2400
    net.Sbase = 1_000_000
    net.Zbase = net.Vbase^2/net.Sbase
    net.bounds.v_lower_mag = 0.95
    net.bounds.v_upper_mag = 1.05
    
    # a radial network has n_edges = n_vertices - 1
    @test net.graph.graph.ne == Graphs.nv(net.graph) - 1  

    net[("650", "rg60")].vreg_pu = dss_voltages["rg60"][1]

    m = build_min_loss_model(net)
    optimize!(m)
    
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

    vsqrd = get_variable_values(:vsqrd, m)
    vs = Dict(k => sqrt.(v) for (k,v) in vsqrd)

    #  I = get_variable_values(:lij, m, p)
    for b in keys(vs)
        try
            @test abs(vs[b][1] - dss_voltages[b][1]) < 0.005
        catch
            println("bus $b failed with difference $(vs[b][1] - dss_voltages[b][1])")
        end
    end

    splitting_bs, subgraph_bs = splitting_busses(net, "650"; max_busses=7)  # 671 and rg60

    mg = split_at_busses(net, splitting_bs, subgraph_bs)

    # test the decomposed solution against openDSS
    builder = Dict(
        v => build_min_loss_model for v in vertices(mg)
    )

    solve_metagraph!(mg, builder, [1e-3, 1e-4, 1e-4]; verbose=false)

    function metagraph_voltages(mg)
        d = Dict()
        for v in vertices(mg)
            merge!(d, get_variable_values(:vsqrd, mg.graph_data[:models][v]))
        end
        return d
    end

    vs = metagraph_voltages(mg)
    for b in keys(vs)
        @test abs(sqrt(vs[b][1]) - dss_voltages[b][1]) < 0.005
    end

end


@testset "basic two-line multiphase" begin
    net = CPF.Network(joinpath("data", "two_line_multi_phase.yaml"))
    net.bounds.i_upper_mag = 15 * net.Sbase / net.Vbase

    m = Model(CSDP.Optimizer)
    set_attribute(m, "printlevel", 0)
    # TODO add pieces of IEEE13 one at a time, starting with two phase lateral with load
    build_bfm!(m, net)

    @objective(m, Min, 
        sum( sum(real.(diag(m[:l][t][i_j]))) for t in 1:net.Ntimesteps, i_j in  edges(net))
    )
    
    optimize!(m)
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]

    # TODO do these vs match the eigen method values? (and should sqrt come last always?)
    # TODO results methods in CommonOPF for multiphase models
    vs = Dict(
        k => sqrt.(abs.(JuMP.value.(w)))
        for (k, w) in m[:w][1]
    )

    Sijs = Dict(
        k => abs.(JuMP.value.(sij))
        for (k, sij) in m[:Sij][1]
    )

    Lijs = Dict(
        k => abs.(JuMP.value.(lij))
        for (k, lij) in m[:l][1]
    )

    # abs(5.6+im*1.2)  # load magnitude per phase

    S0 = abs.(value.(m[:Sj][1][net.substation_bus]))

    # slack bus injection equals flow out on each phase
    for phs in 1:3
        @test S0[phs] ≈ Sijs[("b1", "b2")][phs, phs] atol=1e-4
    end

    # bus b4 has no load so should have approximately same voltage as b2
    for phs in [1, 3]
        @test vs["b2"][phs, phs] ≈ vs["b4"][phs, phs] atol=1e-5
    end

    # no phase 2 at b4
    @test vs["b4"][2, :] == [0,0,0]
    @test vs["b4"][:, 2] == [0,0,0]


    # NLP
    m = Model(Ipopt.Optimizer)
    set_optimizer_attribute(m, "print_level", 0)

    BranchFlowModel.add_bfm_variables(m, net)

    BranchFlowModel.constrain_bfm_nlp(m, net)

    @objective(m, Min, 
        sum( 
            sum(real.(m[:i][t][i_j]) .* real.(m[:i][t][i_j])) 
            + sum(imag.(m[:i][t][i_j]) .* imag.(m[:i][t][i_j])) 
            for t in 1:net.Ntimesteps, i_j in edges(net) 
        )
    )

    optimize!(m)
    
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

    vs = Dict(
        k => abs.(JuMP.value.(v))
        for (k,v) in m[:v][1]
    )

    Sijs = Dict(
        k => abs.(JuMP.value.(sij))
        for (k, sij) in m[:Sij][1]
    )

    # slack bus injection equals flow out on each phase
    for phs in 1:3
        @test S0[phs] ≈ Sijs[("b1", "b2")][phs, phs] atol=1e-4
    end

    # bus b4 has no load so should have approximately same voltage as b2
    for phs in [1, 3]
        @test vs["b2"][phs] ≈ vs["b4"][phs] atol=1e-5
    end

end


@testset "ieee13 unbalanced MultiPhase SDP" begin
    # was testing no load (commented out stuff), now light load but with voltage tol of 2%

    dssfilepath = "data/ieee13/IEEE13_simple_light_load.dss"
    OpenDSS.Text.Command("Redirect $dssfilepath")

    # # set loads to zero
    # kW = 0.0
    # kvar = 0.0
    # for load_name in OpenDSS.Loads.AllNames()
    #     OpenDSS.Loads.Name(load_name)
    #     OpenDSS.CktElement.Enabled(true) # not reached
    #     OpenDSS.Loads.kW(kW)
    #     OpenDSS.Loads.kvar(kvar)
    # end

    OpenDSS.Solution.Solve()

    @test(OpenDSS.Solution.Converged() == true)

    @test(check_opendss_powers() == true)

    dss_voltages = dss_voltages_pu()

    net = CPF.dss_to_Network(dssfilepath)
    
    # for lb in CPF.load_busses(net)
    #     if !(ismissing(net[lb][:Load].kws1))
    #         net[lb][:Load].kws1 = [kW]
    #     end
    #     if !(ismissing(net[lb][:Load].kws2))
    #         net[lb][:Load].kws2 = [kW]
    #     end
    #     if !(ismissing(net[lb][:Load].kws3))
    #         net[lb][:Load].kws3 = [kW]
    #     end
    #     if !(ismissing(net[lb][:Load].kvars1))
    #         net[lb][:Load].kvars1 = [kvar]
    #     end
    #     if !(ismissing(net[lb][:Load].kvars2))
    #         net[lb][:Load].kvars2 = [kvar]
    #     end
    #     if !(ismissing(net[lb][:Load].kvars3))
    #         net[lb][:Load].kvars3 = [kvar]
    #     end
    # end

    net.v0 = 1.02
    net.Vbase = 4160 / sqrt(3)
    net.Sbase = 1e5
    net.Zbase = net.Vbase^2 / net.Sbase
    net.bounds.v_upper_mag = 1.02
    net.bounds.v_lower_mag = 0.95
    net.bounds.s_upper_real =  100
    net.bounds.s_lower_real = 0
    net.bounds.s_upper_imag =  100
    net.bounds.s_lower_imag = 0
    net.bounds.i_upper_mag = 100
    net.bounds.i_lower_mag = 0

    # net[("650", "rg60")].vreg_pu = dss_voltages["rg60"]

    m = Model(CSDP.Optimizer)
    set_attribute(m, "printlevel", 0)

    build_bfm!(m, net)

    @objective(m, Min, 
        sum( sum(real.(diag(m[:l][t][i_j]))) for t in 1:net.Ntimesteps, i_j in edges(net) )
    )

    optimize!(m)
    
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]

    vs = Dict(
        k => diag(sqrt.(abs.(JuMP.value.(w))))
        for (k,w) in m[:w][1]
    )

    @test_nowarn(check_rank_one(m, net))
    
    errors = []
    for b in keys(vs)
        for (i,phsv) in enumerate(filter(v -> v != 0, vs[b]))
            @test abs(phsv - dss_voltages[b][i]) < 2e-2
            # println("$b - $i  $(phsv - dss_voltages[b][i])")
            push!(errors, phsv - dss_voltages[b][i])
        end
    end
    # println(minimum(errors))
    # println(maximum(errors))

end


@testset "SinglePhase network reduction" begin
    # confirm that optimal results do not change

    function make_solve_min_loss_model(net)
        m = build_min_loss_model(net)
        optimize!(m)
        return m
    end

    # 1 validate BFM against OpenDSS
    dssfilepath = "data/singlephase38lines/master.dss"
    net = CPF.dss_to_Network(dssfilepath)
    
    net.bounds.v_upper_mag = 1.05
    net.bounds.v_lower_mag = 0.95

    net.Sbase = 1e5
    net.Vbase = 12.47e3
    net.Zbase = net.Vbase^2 / net.Sbase
    net.v0 = 1.0
    m = make_solve_min_loss_model(net)
    @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

    vsqrd = get_variable_values(:vsqrd, m)
    vs = Dict(k => sqrt.(v) for (k,v) in vsqrd)

    # make the dss solution to compare
    work_dir = pwd()
    OpenDSS.dss("""
        clear
        compile $dssfilepath
        solve
    """)
    cd(work_dir)
    @test(OpenDSS.Solution.Converged() == true)
    dss_voltages = dss_voltages_pu()

    for b in keys(vs)
        @test abs(vs[b][1] - dss_voltages[b][1]) < 0.001
    end

    nvar_original = JuMP.num_variables(m)

    # 2 validate BFM results stay the same after reduction
    nbusses_before = length(busses(net))
#     reduce_tree!(p)
#     m = make_solve_min_loss_model(p)
#     @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]
#     vs_reduced = get_variable_values(:vsqrd, m, p)
#     for b in keys(vs_reduced)
#         @test abs(dss_voltages[b][1] - sqrt(vs_reduced[b][1])) < 0.001
#     end
#     nvar_reduced = JuMP.num_variables(m)
#     @test nvar_original > nvar_reduced  # 225 > 159

#     # 3 split and solve the reduced model, compare v
#     g = BranchFlowModel.make_graph(p.busses, p.edges; directed=true)
#     p_above, p_below = BranchFlowModel.split_inputs(p, "12");
#     # solve above first with sum of p_below loads, then set p_below.v0, solve p_below, set p_above.P/Qload to p_below.substation_bus values
#     # later solve in parallel
#     # p_above.Pload and Qload exist already for bus "12"; we need to add p_below's loads to it
#     init_inputs!([p_below, p_above])  # set p_above loads at p_below.substation_bus
#     # NOTE only one time step so all load vectors have one value
#     @test p_above.Pload["12"][1] == sum( v[1] for v in values(p_below.Pload) )
#     @test p_above.Qload["12"][1] == sum( v[1] for v in values(p_below.Qload) )

#     m_above = make_solve_min_loss_model(p_above)
#     init_vs = Dict(
#         p_below.substation_bus => sqrt(value(m_above[:vsqrd][p_below.substation_bus,1]))
#     )
#     init_inputs!([p_above, p_below]; init_vs=init_vs)
#     @test p_below.v0 == sqrt(value(m_above[:vsqrd][p_below.substation_bus,1]))

#     m_below = make_solve_min_loss_model(p_below)

#     # at this point the v's at bus 12 agree by design
#     # but the loads do not b/c we did not account for losses in first iteration
#     # the plus sign is not a minus because the P/Q values should be equal and _opposite_
#     pdiff = value(m_below[:Pj]["12",1]) + value(m_above[:Pj]["12",1])  # 0.001606
#     qdiff = value(m_below[:Qj]["12",1]) + value(m_above[:Qj]["12",1])  # 0.000526
#     vdiff = 1.0

#     tol = 1e-6
#     while pdiff > tol || qdiff > tol || vdiff > tol
#         p_above.Pload["12"][1] = value(m_below[:Pj]["12",1]) * p_above.Sbase
#         p_above.Qload["12"][1] = value(m_below[:Qj]["12",1]) * p_above.Sbase
#         m_above = make_solve_min_loss_model(p_above)
#         v_above = sqrt(value(m_above[:vsqrd][p_below.substation_bus,1]))
#         vdiff = abs(p_below.v0 - v_above)
#         p_below.v0 = v_above
#         m_below = make_solve_min_loss_model(p_below)
#         pdiff = value(m_below[:Pj]["12",1]) + value(m_above[:Pj]["12",1])  # 3.758466e-8
#         qdiff = value(m_below[:Qj]["12",1]) + value(m_above[:Qj]["12",1])  # 1.231675e-8
#     end

#     vs_decomposed = get_variable_values(:vsqrd, m_above, p_above)
#     merge!(vs_decomposed, get_variable_values(:vsqrd, m_below, p_below))
#     for b in keys(vs_decomposed)
#         @test abs(dss_voltages[b][1] - sqrt(vs_decomposed[b][1])) < 0.001
#     end

#     # split model into 3 models and solve
#     mg = CommonOPF.split_at_busses(p, ["7", "13"])
#     @test mg[1, :p].substation_bus == "0"
#     @test mg[2, :p].substation_bus == "7"
#     @test mg[3, :p].substation_bus == "13"
#     # init_inputs!(mg)  # done in CommonOPF.split_at_busses
#     for v in get_prop(mg, :load_sum_order)
#         set_prop!(mg, v, :m, make_solve_min_loss_model(mg[v, :p]))
#     end
#     set_indexing_prop!(mg, :m)
#     # every model solved once using load approximations
#     # now set better load approximations and connect voltages
#     set_inputs!(mg)
#     @test mg[2,:p].v0[1] ≈ sqrt(value(mg[1,:m][:vsqrd][mg[2,:p].substation_bus,1]))
#     @test mg[3,:p].v0[1] ≈ sqrt(value(mg[2,:m][:vsqrd][mg[3,:p].substation_bus,1]))
#     # run another solve
#     for v in get_prop(mg, :load_sum_order)
#         set_prop!(mg, v, :m, make_solve_min_loss_model(mg[v, :p]))
#     end
#     pdiffs1, qdiffs1, vdiffs1 = get_diffs(mg)
#     # another round to compare:
#     set_inputs!(mg)
#     for v in get_prop(mg, :load_sum_order)
#         set_prop!(mg, v, :m, make_solve_min_loss_model(mg[v, :p]))
#     end
#     pdiffs2, qdiffs2, vdiffs2 = get_diffs(mg)
#     @test sum(pdiffs2) < sum(pdiffs1) 
#     @test sum(qdiffs2) < sum(qdiffs1) 
#     @test sum(vdiffs2) < sum(vdiffs1) 

#     vs_decomposed = get_variable_values(:vsqrd, mg[1, :m], mg[1, :p])
#     merge!(vs_decomposed, get_variable_values(:vsqrd, mg[2, :m], mg[2, :p]))
#     merge!(vs_decomposed, get_variable_values(:vsqrd, mg[3, :m], mg[3, :p]))
#     for b in keys(vs_decomposed)
#         @test abs(dss_voltages[b][1] - sqrt(vs_decomposed[b][1])) < 0.001
#     end
end


# @testset "Single phase regulators at network splits" begin
#     Sbase = 1e6
#     Vbase = 12.47e3
#     # use OpenDSS to replace a line 20-21 with a transformer/regulator and test solution convergence
#     p = Inputs(
#         joinpath("data", "singlephase38lines", "master_extra_trfx.dss"), 
#         "0";
#         Sbase=Sbase, 
#         Vbase=Vbase, 
#         v0 = 1.00,
#        .bounds.v_upper_mag = 1.05,
#        .bounds.v_lower_mag = 0.95,
#         relaxed = false,
#     );
#     @test ("20","21") in keys(p.regulators)
#     @test reg_busses(p) == ["21"]

#     dss("clear")
#     dss("Redirect data/singlephase38lines/master_extra_trfx.dss")
#     dss("Solve")
#     @test(OpenDSS.Solution.Converged() == true)
#     dss_voltages = dss_voltages_pu()

#     p.regulators[("20","21")][:turn_ratio] = dss_voltages["21"][1] / dss_voltages["20"][1]
#     @test turn_ratio(p, "21") == dss_voltages["21"][1] / dss_voltages["20"][1]

#     # split at two busses including the regulator and compare against openDSS
#     mg = CommonOPF.split_at_busses(p, ["14","21"])
#     builder = Dict(
#         v => build_min_loss_model for v in vertices(mg)
#     )
#     solve_metagraph!(mg, builder, [1e-5, 1e-5, 1e-5]; verbose=false)
#     vs = metagraph_voltages(mg)
#     for b in keys(vs)
#         @test abs(sqrt(vs[b][1]) - dss_voltages[b][1]) < 0.001
#     end

#     # split only at regulator, add vreg, and test the vdiff is zero
#     # (because using vreg makes the regulator voltage equal to the setting)
#     p.regulators[("20","21")][:vreg] = 1.02
#     mg = CommonOPF.split_at_busses(p, ["21"])

#     builder = Dict(
#         v => build_min_loss_model for v in vertices(mg)
#     )
#     solve_metagraph!(mg, builder, [1e-3, 1e-3, 1e-3]; verbose=false)
#     pdiffs, qdiffs, vdiffs = get_diffs(mg)
#     @test maximum(vdiffs) ≈ 0  # because that is how it is defined across a regulator


# end

end  # all tests
