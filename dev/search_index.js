var documenterSearchIndex = {"docs":
[{"location":"methods/#Methods","page":"Methods","title":"Methods","text":"","category":"section"},{"location":"methods/","page":"Methods","title":"Methods","text":"Some various methods used in BranchFlowModel.jl:","category":"page"},{"location":"methods/","page":"Methods","title":"Methods","text":"warning: Warning\nThis list of exported methods may not be up to date and there are missing doc strings. Contributions are welcome via fork and pull request.","category":"page"},{"location":"methods/","page":"Methods","title":"Methods","text":"Inputs\nsinglephase38linesInputs\ndsstxt_to_sparse_array \nparse_dss\nbuild_model!\nadd_variables\nconstrain_power_balance\nconstrain_substation_voltage\nconstrain_KVL\nconstrain_loads\nconstrain_bounds\ncheck_rank_one\nget_bus_values \nget_edge_values \nget_ijlinecode\nget_ij_idx\n# recover_voltage_current  # TODO validate this method\ni_to_j \nj_to_k \nrij \nxij\nzij\ncheck_soc_inequalities\ncheck_connected_graph\nget_load_bal_shadow_prices\nvoltage_values_by_time_bus\ncurrent_values_by_time_edge\nline_flow_values_by_time_edge\nreduce_tree!\ntrim_tree!\nmake_graph\nleaf_busses\ninit_inputs!\nset_inputs!\nsplit_at_busses\nget_diffs\nall_outneighbors\nall_inneighbors\nbusses_from_deepest_to_source\nvertices_from_deepest_to_source\nsplitting_busses\nconnect_subgraphs_at_busses\nsplit_inputs\nbuild_metagraph\nsolve_metagraph!\nmetagraph_voltages\ncheck_unique_solution_conditions","category":"page"},{"location":"methods/#BranchFlowModel.singlephase38linesInputs","page":"Methods","title":"BranchFlowModel.singlephase38linesInputs","text":"singlephase38linesInputs(;\n    Pload=Dict{String, AbstractArray{Real, 1}}(), \n    Qload=Dict{String, AbstractArray{Real, 1}}(), \n    T=24,\n    loadnodes = [\"3\", \"5\", \"36\", \"9\", \"10\", \"11\", \"12\", \"13\", \"15\", \"17\", \"18\", \"19\", \"22\", \"25\", \n                \"27\", \"28\", \"30\", \"31\", \"32\", \"33\", \"34\", \"35\"],\n    Sbase = 1e6,\n    Vbase = 12.5e3,\n    v0=1.0,\n    v_uplim = 1.05,\n    v_lolim = 0.95,\n)\n\nConvenience function for creating a single phase network with 38 lines and nodes.  Taken from: Andrianesis et al. 2019 \"Locational Marginal Value of Distributed Energy Resources as Non-Wires Alternatives\"\n\nNOTE that Inputs is a mutable struct (s.t. loads can be added later).\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.dsstxt_to_sparse_array","page":"Methods","title":"BranchFlowModel.dsstxt_to_sparse_array","text":"dsstxt_to_sparse_array(fp::String, first_data_row::Int = 5)\n\nconvert a SystemY.txt file from OpenDSS to a julia matrix. assumes that Y is symmetric.\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.add_variables","page":"Methods","title":"BranchFlowModel.add_variables","text":"add_variables(m, p::Inputs{MultiPhase})\n\nCreate complex variables:\n\nm[:w] are 3x3 Hermitian matrices of voltage squared (V*V^T)\nm[:l] are 3x3 Hermitian matrices of current squared (I*I^T)\nm[:Sj] are 3x1 matrices of net power injections (at bus j)\nm[:Sij] are 3x3 Complex matrices of line flow powers (from i to j)\n\nThe positive semi-definite constraints are also defined and stored as\n\nm[:H][t][j] where t is time step and j is the bus name.\n\nAll of the variable containers have typeof Dict{Int, Dict{String, AbstractVecOrMat}}`.\n\nThe first index is time step (integer)\nThe second index is bus or line (string)\nand finally a matrix of complex variables\n\nSome examples of using variables:\n\n\nvalue.(m[:Sj][1][\"671\"])\n\nvalue(variable_by_name(m, \"real(Sj_1_671_1)\"))\n\nvalue.(m[:w][1][\"671\"])\n\nvalue.(m[:Sj][1][p.substation_bus]) \n\nfor b in keys(p.Pload)\n    println(b, \"  \", value.(m[:Sj][1][b]))\nend\n\nfix(variable_by_name(m, \"real(Sj_1_645_3)\"), 0.0, force=true)\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.constrain_power_balance","page":"Methods","title":"BranchFlowModel.constrain_power_balance","text":"function constrain_power_balance(m, p::Inputs)\n\nPij in - losses == sum of line flows out + net injection NOTE: using sum over Pij for future expansion to mesh grids i -> j -> k\n\n\n\n\n\nfunction constrain_power_balance(m, p::Inputs{MultiPhase})\n\nSij in - losses == sum of line flows out + net injection NOTE: using sum over Pij for future expansion to mesh grids i -> j -> k\n\nAll of the power balance constraints are stored in m[:loadbalcons] with the bus name (string) as the first index. For example m[:loadbalcons][\"busname\"] will give the constrain container from JuMP for all time steps.\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.constrain_KVL","page":"Methods","title":"BranchFlowModel.constrain_KVL","text":"constrain_KVL(m, p::Inputs{MultiPhase})\n\nAdd the voltage drop definintions between busses.\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.constrain_loads","page":"Methods","title":"BranchFlowModel.constrain_loads","text":"constrain_loads(m, p::Inputs)\n\nset net injections Pj/Qj to negative of Inputs.Pload/Qload, which are normalized by Sbase when creating Inputs\nkeys of P/Qload must match Inputs.busses. Any missing keys have load set to zero.\nInputs.substation_bus is unconstrained, slack bus\n\n\n\n\n\nconstrain_loads(m, p::Inputs{MultiPhase})\n\nset loads to negative of Inputs.Pload and Inputs.Qload,    which are normalized by Sbase when creating Inputs.\nkeys of Pload and Qload must match Inputs.busses. Any missing keys have load set to zero.\nInputs.substation_bus is unconstrained, slack bus\n\nEach of the power injection constraints are stored in the model under m[:injectioncons]. To acces the constraints you can:\n\nm[:injectioncons][\"busname\"][t,phs]\n\nwhere t is the integer time step and phs is the integer phase.\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.check_rank_one","page":"Methods","title":"BranchFlowModel.check_rank_one","text":"check_rank_one(m::JuMP.AbstractModel, p::Inputs{BranchFlowModel.MultiPhase}, tol=1e-3)\n\nCheck the rank of the m[:H] matrices from the PSD cone constraints. Warnings express any values with rank greater than one.\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.i_to_j","page":"Methods","title":"BranchFlowModel.i_to_j","text":"function i_to_j(j::AbstractString, p::Inputs)\n\nfind all busses upstream of bus j\n\nnote: Note\nIn a radial network this function should return an Array with length of 1.\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.j_to_k","page":"Methods","title":"BranchFlowModel.j_to_k","text":"function j_to_k(j::AbstractString, p::Inputs)\n\nfind all busses downstream of bus j\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.rij","page":"Methods","title":"BranchFlowModel.rij","text":"rij(i::AbstractString, j::AbstractString, p::Inputs{SinglePhase})\n\nThe per-unit resistance of line i->j\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.xij","page":"Methods","title":"BranchFlowModel.xij","text":"xij(i::AbstractString, j::AbstractString, p::Inputs{SinglePhase})\n\nThe per-unit reacttance of line i->j\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.check_soc_inequalities","page":"Methods","title":"BranchFlowModel.check_soc_inequalities","text":"check_soc_inequalities(m::JuMP.AbstractModel, p::Inputs)\n\ncreate and return a vector of the gaps in the second order cone constraints\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.check_connected_graph","page":"Methods","title":"BranchFlowModel.check_connected_graph","text":"check_connected_graph(p::Inputs{BranchFlowModel.SinglePhase})\n\nreturn true if only one connected graph; false otherwise\n\nthis is a good check to do before attempting to solve a model because if there is more than one sub-graph then it is likely the model will be infeasible\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.get_load_bal_shadow_prices","page":"Methods","title":"BranchFlowModel.get_load_bal_shadow_prices","text":"get_load_bal_shadow_prices(m::JuMP.AbstractModel, p::Inputs)\n\ncreate and return a dict indexed by bus and time for shadow prices     (just real prices for now)\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.reduce_tree!","page":"Methods","title":"BranchFlowModel.reduce_tree!","text":"reduce_tree!(p::Inputs{SinglePhase})\n\ncombine any line sets with intermediate busses that have indegree == outdegree == 1 and is not a load bus into a single line\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.make_graph","page":"Methods","title":"BranchFlowModel.make_graph","text":"make_graph(busses::AbstractVector{String}, edges::AbstractVector{Tuple})\n\nreturn SimpleDiGraph, Dict, Dict  with the dicts for bus => int and int => bus (because Graphs.jl only works with integer nodes)\n\njulia> g[\"13\", :bus]\n10\n\njulia> get_prop(g, :bus_int_map)[\"13\"]\n10\n\njulia> g[13, :bus]\n\"24\"\n\njulia> get_prop(g, :int_bus_map)[13]\n\"24\"\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.leaf_busses","page":"Methods","title":"BranchFlowModel.leaf_busses","text":"leaf_busses(p::Inputs)\n\nreturns Vector{String} containing all of the leaf busses in p.busses\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.init_inputs!","page":"Methods","title":"BranchFlowModel.init_inputs!","text":"init_inputs!(ps::Vector{Inputs{BranchFlowModel.SinglePhase}}; init_vs::Dict = Dict())\n\nSet the load on the upstream leaf noades equal to the sum of all the loads in the downstream inputs. It is important that the order of ps is from leaf branches to trunk branches so that the sums of loads take into account all downstream sub-trees.\n\nif init_vs is provided, the p.v0 is set for the Input with its p.substation_bus equal      to the key in init_vs\n\ninit_vs = Dict(\n    \"sub_bus_1\" => 0.98\n)\n\nfor p in ps\n    if p.substation_bus in keys(init_vs)\n        p.v0 = init_vs[p.substation_bus]\n    end\nend\n\n\n\n\n\ninit_inputs!(mg::MetaDiGraph; init_vs::Dict = Dict())\n\nUse the :load_sum_order in mg to init_inputs! in the correct order, i.e. set the loads at the leaf - substation connections as sums of all the loads (and the voltages at substations)\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.set_inputs!","page":"Methods","title":"BranchFlowModel.set_inputs!","text":"set_inputs!(mg::MetaDiGraph; α::Float64=0.0)\n\nSet the shared values in each subgraph / vertex of mg:\n\nset the current vertex's v0 to its inneighbor's voltage\nset the current vertex P/Qload to the outneighbors' substation_bus loads\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.split_at_busses","page":"Methods","title":"BranchFlowModel.split_at_busses","text":"split_at_busses(p::Inputs{BranchFlowModel.SinglePhase}, at_busses::Vector{String})\n\nSplit Inputs using the at_busses\n\nreturns MetaDiGraph with vertex properties :p containing Input for the sub-graphs. For example mg[2, :p] is the Input at the second vertex of the graph created by splitting  the network via the at_busses.\n\n\n\n\n\nsplit_at_busses(p::Inputs{BranchFlowModel.SinglePhase}, at_busses::Vector{String}, with_busses::Vector{Vector{String}})\n\nSplit up p using the at_busses as each new substation_bus and containing the corresponding with_busses. The at_busses and with_busses can be determined using splitting_busses.\n\nNOTE: this variation of spltatbusses allows for more than two splits at the same bus; whereas the other implementation of splitatbusses only splits the network into two parts for everything above and everything below a splitting bus.\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.get_diffs","page":"Methods","title":"BranchFlowModel.get_diffs","text":"get_diffs(mg::MetaDiGraph)\n\nUses the JuMP Models stored in mg[:m] to calculate the difference between Pj, Qj, and |v| at every leaf/substation connection. \n\nreturns three Float64[]\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.busses_from_deepest_to_source","page":"Methods","title":"BranchFlowModel.busses_from_deepest_to_source","text":"busses_from_deepest_to_source(g::MetaDiGraph, source::String)\n\nreturn the busses and their integer depths in order from deepest from shallowest\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.splitting_busses","page":"Methods","title":"BranchFlowModel.splitting_busses","text":"splitting_busses(p::Inputs{BranchFlowModel.SinglePhase}, source::String; threshold::Int64=10)\n\nDetermine the busses to split a tree graph on by searching upward from the deepest leafs first and gathering the nearest busses until threshold is met for each subgraph.\n\nReturns a Vector{String} for the bus names and Vector{Vector{String}} for the corresponding busses within each sub-graph.\n\nnote: Note\nIt is not enough to have only the splitting busses to obey the max_busses limit because one must also know which sub branches to take from each splitting bus. In other words, we also need all the busses within each subgraph to split properly. For example, if a splitting bus has two sub branches then obeying the max_busses limit can require only including one sub branch out of the splitting bus. To know which branch to take we can use the other busses in the sub graph.\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.connect_subgraphs_at_busses","page":"Methods","title":"BranchFlowModel.connect_subgraphs_at_busses","text":"connect_subgraphs_at_busses(p::Inputs{BranchFlowModel.SinglePhase}, at_busses::Vector{String}, subgraphs::Vector{Vector})\n\nThe splitting_busses algorithm does not include over laps in subgraphs. But, we want overlaps at the splitting busses for solving the decomposed branch flow model. So here we add the overlapping splitting busses to each sub graph.\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.split_inputs","page":"Methods","title":"BranchFlowModel.split_inputs","text":"split_inputs(p::Inputs{BranchFlowModel.SinglePhase}, bus::String, g::SimpleDiGraph)\n\nSplit inputs into one graph for everything above bus and one graph for everything     below bus.\n\n\n\n\n\nsplit_inputs(p::Inputs{BranchFlowModel.SinglePhase}, bus::String, out_buses::Vector{String})\n\nSplit p into p_above and p_below where p_below has only out_buses and p_above has union( [bus], setdiff(p.busses, out_buses) ).\n\nNote that out_buses must contain bus\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.build_metagraph","page":"Methods","title":"BranchFlowModel.build_metagraph","text":"build_metagraph(p::Inputs{BranchFlowModel.SinglePhase}, source::String; max_busses::Int64=10)\n\nreturn MetaDiGraph with :p property set for every vertex by splitting the Inputs via splitting_busses\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.solve_metagraph!","page":"Methods","title":"BranchFlowModel.solve_metagraph!","text":"solve_metagraph!(mg::MetaDiGraph, builder::Function, tol::T; α::T=0.5, verbose=false) where T <: Real\n\nGiven a MetaDiGraph and a JuMP Model builder method iteratively solve the models until the tol is  met for the differences provided by BranchFlowModel.get_diffs. \n\nThe builder must accept only one argument of type BranchFlowModel.AbstractInputs that returns  a JuMP.AbstractModel. Each model returned from the builder is stored as an :m property in  each vertex of mg.\n\nnote: Note\ntol is compared to the maximum absolute value of all the p, q, and v differences.\n\n\n\n\n\nsolve_metagraph!(mg::MetaDiGraph, builder::Function, tols::Vector{T}; α::T=0.5, verbose=false) where T <: Real\n\nGiven a MetaDiGraph and a JuMP Model builder method iteratively solve the models until the tols are  met for the differences provided by BranchFlowModel.get_diffs. \n\nThe builder must accept only one argument of type BranchFlowModel.AbstractInputs that returns  a JuMP.AbstractModel. Each model returned from the builder is stored as an :m property in  each vertex of mg.\n\nnote: Note\nThe tols should have a length of three. The first value is compared to the maximum absolute difference in Pj, the second for Qj, and the third for |v|. All differences are calculated at the leaf/substation connections.\n\n\n\n\n\nsolve_metagraph!(mg::MetaDiGraph, builder::Dict{Int64, Function}, tols::Vector{T}; α::T=0.5, verbose=false) where T <: Real\n\nGiven a MetaDiGraph and a JuMP Model builder method iteratively solve the models until the tols are  met for the differences provided by BranchFlowModel.get_diffs.  The builder dict is used to build each model for the corresponding vertex key.\n\nEach function in the builder dict must accept only one argument of type BranchFlowModel.AbstractInputs that returns  a JuMP.AbstractModel. Each model returned from the builder function is stored as an :m property in  each vertex of mg.\n\nnote: Note\nThe tols should have a length of three. The first value is compared to the maximum absolute difference in Pj, the second for Qj, and the third for |v|. All differences are calculated at the leaf/substation connections.\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.check_unique_solution_conditions","page":"Methods","title":"BranchFlowModel.check_unique_solution_conditions","text":"check_unique_solution_conditions(p::Inputs)\n\nreport the maximum per-unit immpedance and load values. See Chiang and Baran 2013: A load flow solution with feasible voltage magnitude always exists and is unique when\n\nV0 ≈ 1\nloss values < 1\nrpu, xpu << 1\n\n\n\n\n\n","category":"function"},{"location":"math/#Single-Phase-BranchFlowModel","page":"Math","title":"Single Phase BranchFlowModel","text":"","category":"section"},{"location":"math/","page":"Math","title":"Math","text":"From [1]","category":"page"},{"location":"math/","page":"Math","title":"Math","text":"Notation:","category":"page"},{"location":"math/","page":"Math","title":"Math","text":"P_ij real power flow from node i to node j\np_j real power injection on node j\n`\\\\mathcal{N}^+ set of all nodes in network except the source\nw_j voltage magnitude squared on node j\nell_ij current magnitude squared  from node i to node j","category":"page"},{"location":"math/","page":"Math","title":"Math","text":"beginaligned\nP_ij - r_ij ell_ij + p_j = sum_kjrightarrow k P_jk  forall j in mathcalN^+ \nQ_ij - x_ij ell_ij + q_j = sum_kjrightarrow k Q_jk  forall j in mathcalN^+ \nw_j = w_i - 2 r_ij P_ij - 2 x_ij Q_ij + (r_ij^2 + x_ij^2) ell_ij  forall j in mathcalN^+ \nw_i ell_ij = P_ij^2 + Q_ij^2 forall (ij) in mathcalE \n(v_jmin)^2 le w_j le (v_jmax)^2  forall j in mathcalN^+ \nendaligned","category":"page"},{"location":"math/#Three-Phase-BranchFlowModel","page":"Math","title":"Three Phase BranchFlowModel","text":"","category":"section"},{"location":"math/","page":"Math","title":"Math","text":"TODO","category":"page"},{"location":"math/#References","page":"Math","title":"References","text":"","category":"section"},{"location":"math/#[1]","page":"Math","title":"[1]","text":"","category":"section"},{"location":"math/","page":"Math","title":"Math","text":"Baran, Mesut E., and Felix F. Wu. \"Optimal capacitor placement on radial distribution systems.\" IEEE Transactions on power Delivery 4.1 (1989): 725-734. Chicago\t","category":"page"},{"location":"#BranchFlowModel.jl","page":"User Documentation","title":"BranchFlowModel.jl","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"BranchFlowModel builds the branch flow constraints using JuMP.  The intent of this package is to allow users to build mathematical programs that include BranchFlowModel constraints. No objective is added to the JuMP model in this package and so solving any problem defined by the constraints built by BranchFlowModel.jl is a feasibility problem. Dictionaries of constraints are provided so that one can delete and/or modify the base constraints to fit their problem.","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"warning: Warning\nThis package is under development. Contributions are welcome via fork and pull request.","category":"page"},{"location":"#Inputs","page":"User Documentation","title":"Inputs","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"There are two methods for creating Inputs:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Using openDSS files\nProviding the network topology","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Inputs(::String, ::String)\nInputs(::AbstractVector{<:Tuple}, ::AbstractVector{<:AbstractString}, ::AbstractVector{<:Real}, ::AbstractVector{<:AbstractVector}, ::String)","category":"page"},{"location":"#BranchFlowModel.Inputs-Tuple{String, String}","page":"User Documentation","title":"BranchFlowModel.Inputs","text":"Inputs(\n    dssfilepath::String, \n    substation_bus::String;\n    Pload::AbstractDict=Dict(), \n    Qload::AbstractDict=Dict(), \n    Sbase=1, \n    Vbase=1, \n    v0, \n    v_lolim=0.95, \n    v_uplim=1.05,\n    Ntimesteps=1, \n    P_up_bound=1e4,\n    Q_up_bound=1e4,\n    P_lo_bound=-1e4,\n    Q_lo_bound=-1e4,\n    relaxed=true,\n)\n\nInputs constructor that parses a openDSS file for the network. If Pload and Qload are not provided then the loads are also parsed from the openDSS file.\n\n\n\n\n\n","category":"method"},{"location":"#BranchFlowModel.Inputs-Tuple{AbstractVector{<:Tuple}, AbstractVector{<:AbstractString}, AbstractVector{<:Real}, AbstractVector{<:AbstractVector}, String}","page":"User Documentation","title":"BranchFlowModel.Inputs","text":"Inputs(\n    edges::Array{Tuple}, \n    linecodes::Array{String}, \n    linelengths::Array{Float64}, \n    phases::Vector{Vector},\n    substation_bus::String;\n    Pload, \n    Qload, \n    Sbase=1, \n    Vbase=1, \n    Zdict, \n    v0, \n    v_lolim=0.95, \n    v_uplim=1.05,\n    Ntimesteps=1, \n    P_up_bound=1e4,\n    Q_up_bound=1e4,\n    P_lo_bound=-1e4,\n    Q_lo_bound=-1e4,\n    Isquared_up_bounds=Dict{String, Float64}(),\n    relaxed=true\n)\n\nLowest level Inputs constructor (the only one that returns the Inputs struct). \n\nnote: Note\nThe real and reactive loads provided are normalized using Sbase.\n\n\n\n\n\n","category":"method"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Both of the Inputs functions return a mutable Inputs struct:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Inputs","category":"page"},{"location":"#BranchFlowModel.Inputs","page":"User Documentation","title":"BranchFlowModel.Inputs","text":"mutable struct Inputs{T<:Phases} <: AbstractInputs\n    edges::Array{Tuple, 1}\n    linecodes::Array{String, 1}\n    linelengths::Array{Float64, 1}\n    busses::Array{String}\n    phases::Vector{Vector}\n    substation_bus::String\n    Pload::Dict{String, Any}\n    Qload::Dict{String, Any}\n    Sbase::Real\n    Vbase::Real\n    Ibase::Real\n    Zdict::Dict{String, Dict{String, Any}}\n    v0::Real\n    v_lolim::Real\n    v_uplim::Real\n    Zbase::Real\n    Ntimesteps::Int\n    pf::Float64\n    Nnodes::Int\n    P_up_bound::Float64\n    Q_up_bound::Float64\n    P_lo_bound::Float64\n    Q_lo_bound::Float64\n    Isquared_up_bounds::Dict{String, <:Real}\n    phases_into_bus::Dict{String, Vector{Int}}\nend\n\nInputs\n\nedges Vector{Tuple} e.g. [(\"0\", \"1\"), (\"1\", \"2\")]\nlinecodes vector of string keys for the Zdict (impedance values for lines). When using an OpenDSS model a linecode is the name in New linecode.name\nlinelengths vector of floats to scale impedance values\nbusses vector of bus names\nphases vector of vectors with ints for the line phases (e.g. [[1,2,3], [1,3], ...])\nPload dict with busses for keys and uncontrolled real power loads (positive is load) by phase and time\nQload dict with busses for keys and uncontrolled reactive power loads (positive is load) by phase and time\nSbase base apparent power for network, typ. feeder capacity. Used to normalize all powers in model\nVbase base voltage for network, used to determine Zbase = Vbase^2 / Sbase\nIbase = Sbase / (Vbase * sqrt(3))\nZdict dict with linecodes for keys and subdicts with \"xmatrix\" and \"zmatrix\" keys with per unit length values. Values are divided by Zbase and multiplied by linelength in mathematical model.\nv0 slack bus reference voltage\n\nTODO Zdict example\n\nTODO test against simple model to make sure scaling is done right\n\nnote: Note\nThe edges, linecodes, phases, edge_keys, and linelengths are in mutual order (i.e. the i-th value in each list corresponds to the same line)\n\n\n\n\n\n","category":"type"},{"location":"#Building-a-Model","page":"User Documentation","title":"Building a Model","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"The build_ldf! function takes a JuMP.Model and Inputs struct as its two arguments and adds the variables and constraints:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"build_model!","category":"page"},{"location":"#BranchFlowModel.build_model!","page":"User Documentation","title":"BranchFlowModel.build_model!","text":"build_model!(m::JuMP.AbstractModel, p::Inputs)\n\nAdd variables and constraints to m using the values in p. Calls the following functions:\n\nadd_variables(m, p)\nconstrain_power_balance(m, p)\nconstrain_substation_voltage(m, p)\nconstrain_KVL(m, p)\nconstrain_cone(m, p)\nconstrain_loads(m, p)\n\n\n\n\n\nbuild_model!(m::JuMP.AbstractModel, p::Inputs{MultiPhase})\n\nAdd variables and constraints to m using the values in p. Calls the following functions:\n\nadd_variables(m, p)\nconstrain_power_balance(m, p)\nconstrain_substation_voltage(m, p)\nconstrain_KVL(m, p)\nconstrain_loads(m, p)\n\n\n\n\n\n","category":"function"},{"location":"#Variables","page":"User Documentation","title":"Variables","text":"","category":"section"},{"location":"#Single-Phase-Model","page":"User Documentation","title":"Single Phase Model","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Let m be the JuMP.Model provided by the user, then the variables can be accessed via:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"m[:vsqrd] voltage magnitude squared, indexed on busses, (phases), time\nm[:Pj], m[:Qj] net real, reactive power injection, indexed on busses, (phases), time\nm[:Pij], m[:Qij] net real, reactive line flow, indexed on edges, (phases), time","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Note that the m[:Pj], m[:Qj] are not truly variables since they are defined by the loads (unless you modify the model to make them decisions). After a model has been solved using JuMP.optimize! variable values can be extracted with JuMP.value. For more see Getting started with JuMP.","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"note: Note\nSingle phase models do not have a phase index","category":"page"},{"location":"#MultiPhase-Model","page":"User Documentation","title":"MultiPhase Model","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"The definition of the muliphase variables is done in model_multi_phase.jl as follows:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"# voltage squared is Hermitian\nm[:w] = Dict{Int64, S}()\n# current squared is Hermitian\nm[:l] = Dict{Int64, S}()\n# complex line powers (at the sending end)\nm[:Sij] = Dict{Int64, S}()\n# complex net powers injections \nm[:Sj] = Dict{Int64, S}()\n# Hermitian PSD matrices\nm[:H] = Dict{Int64, S}()","category":"page"},{"location":"#Accessing-and-Modifying-Constraints","page":"User Documentation","title":"Accessing and Modifying Constraints","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Let the JuMP.Model provided by the user be called m.  Some constraints are stored in the model dict as anonymous constraints with symbol keys.","category":"page"},{"location":"#Power-Injections","page":"User Documentation","title":"Power Injections","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"BranchFlowModel.jl uses the convention that power injections are positive (and loads are negative). If no load is provided for a given bus (and phase) then the real and reactive power injections at that bus (and phase) are set to zero with an equality constraint.","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"All power injection constraints are stored in m[:injection_equalities]. The constraints are indexed in the following order:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"by bus name (string), as provided in Inputs.busses;\nby \"p\" or \"q\" for real and reactive power respectively;\nby phase number (integer); and\nby time (integer).","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"For example, m[:injection_equalities][\"680\"][\"p\"][2][1] contains the constraint reference for the power injection equality constraint for bus \"680\", real power, in time step 2, on phase 1.","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"If one wished to replace any constraint one must first delete the constraint using the delete function. For example:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"delete(m, m[:cons][:injection_equalities][\"680\"][\"p\"][1])","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Note that the time index was not provided in the delete command in this example, which implies that the equality constraints for all time steps were deleted. One can also delete individual time step constraints by providing the time index.","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"The deleted constraints can then be replaced with a new set of constraints. For example:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"m[:cons][:injection_equalities][\"680\"][\"p\"][1] = @constraint(m, [t in 1:p.Ntimesteps],\n    m[:Pj][\"680\",1,t] == -1e3 / p.Sbase\n)","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"where p is short for \"parameters\" and is the Inputs struct for the problem of interest. Note that it is not necessary to store the new constraints in the m[:cons][:injection_equalities].","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"See the JuMP documentation for more on deleting constraints.","category":"page"}]
}
