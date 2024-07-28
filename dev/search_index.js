var documenterSearchIndex = {"docs":
[{"location":"decomposition/","page":"Decomposition","title":"Decomposition","text":"Methods to decompose and solve the SinglePhase Branch Flow Model are provided based on the work in  [2]. These methods are most advantageous when solving the non-linear (unrelaxed) power flow equations and are only valid(?) in radial networks.","category":"page"},{"location":"decomposition/","page":"Decomposition","title":"Decomposition","text":"init_inputs!\nset_inputs!\nget_diffs\nsplit_inputs","category":"page"},{"location":"decomposition/#BranchFlowModel.set_inputs!","page":"Decomposition","title":"BranchFlowModel.set_inputs!","text":"set_inputs!(mg::MetaGraphsNext.MetaGraph; α::Float64=0.0)\n\nSet the shared values in each subgraph / vertex of mg:\n\nset the current vertex's v0 to its inneighbor's voltage\nset the current vertex P/Qload to the outneighbors' substation_bus loads\n\n\n\n\n\n","category":"function"},{"location":"decomposition/#BranchFlowModel.get_diffs","page":"Decomposition","title":"BranchFlowModel.get_diffs","text":"get_diffs(mg::MetaGraphsNext.MetaGraph)\n\nUses the JuMP Models stored in mg.graph_data[:models] to calculate the difference between power injections/loads, and |v| at every leaf/substation connection. \n\nreturns three Float64[]\n\n\n\n\n\n","category":"function"},{"location":"decomposition/#[2]","page":"Decomposition","title":"[2]","text":"","category":"section"},{"location":"decomposition/","page":"Decomposition","title":"Decomposition","text":"Sadnan, Rabayet, and Anamika Dubey. \"Distributed optimization using reduced network equivalents for radial power distribution systems.\" IEEE Transactions on Power Systems 36.4 (2021): 3645-3656.","category":"page"},{"location":"methods/#Methods","page":"Methods","title":"Methods","text":"","category":"section"},{"location":"methods/","page":"Methods","title":"Methods","text":"Some various methods used in BranchFlowModel.jl:","category":"page"},{"location":"methods/","page":"Methods","title":"Methods","text":"warning: Warning\nThis list of exported methods may not be up to date and there are missing doc strings. Contributions are welcome via fork and pull request.","category":"page"},{"location":"methods/","page":"Methods","title":"Methods","text":"dsstxt_to_sparse_array \ndss_files_to_dict\nbuild_model!(m::JuMP.AbstractModel, net::Network{SinglePhase}; relaxed::Bool=true)build_model!(m::JuMP.AbstractModel, net::Network{MultiPhase}; PSD::Bool=true)add_variables\nconstrain_power_balance\nconstrain_substation_voltage\nconstrain_KVL\nconstrain_bounds\ncheck_rank_one\nget_bus_values \nget_edge_values \ncheck_soc_inequalities\nget_load_bal_shadow_prices\ncurrent_values_by_time_edge\nline_flow_values_by_time_edge\nreduce_tree!\ntrim_tree!\nmake_graph\nleaf_busses\nset_inputs!\nget_diffs\nsolve_metagraph!\nmetagraph_voltages\ncheck_unique_solution_conditions\ncheck_statuses\nreg_busses\nremove_bus!","category":"page"},{"location":"methods/#BranchFlowModel.constrain_power_balance","page":"Methods","title":"BranchFlowModel.constrain_power_balance","text":"function constrain_power_balance(m, net::Network)\n\nDefine the m[:loadbalcons][bus] ∀ bus ∈ busses(net) as a Dict of constraints. The keys are \"p\" and \"q\" for real and reactive power balance respectively. The values are the JuMP constraints.\n\n∑ Pij in - losses + net injection - ∑ Pjk out = 0\n\nThe net injection are user defined loads. If one wishes to make the net injection a decision variable then delete the constraint and redefine the constraint with your decision variable.\n\nNOTE: using sum over Pij for future expansion to mesh grids and the convention: i -> j -> k\n\n\n\n\n\nfunction constrain_power_balance(m, net::Network{MultiPhase})\n\nSij in - losses == sum of line flows out + net injection NOTE: using sum over Pij for future expansion to mesh grids i -> j -> k\n\nAll of the power balance constraints are stored in m[:loadbalcons] with the bus name (string) as the first index. For example m[:loadbalcons][\"busname\"] will give the constrain container from JuMP for all time steps.\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.constrain_KVL","page":"Methods","title":"BranchFlowModel.constrain_KVL","text":"constrain_KVL(m, net::Network{MultiPhase})\n\nAdd the voltage drop definitions between busses.\n\nw_j = w_i - S_ij Z^star - Z S_ij^star + Z L_ij Z^star\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.check_rank_one","page":"Methods","title":"BranchFlowModel.check_rank_one","text":"check_rank_one(m::JuMP.AbstractModel, net::Network, tol=1e-3)\n\nCheck the rank of the m[:H] matrices from the PSD cone constraints. Warnings express any values with rank greater than one.\n\n\n\n\n\n","category":"function"},{"location":"methods/#CommonOPF.reduce_tree!","page":"Methods","title":"CommonOPF.reduce_tree!","text":"reduce_tree!(net::Network{SinglePhase})\n\ncombine any line sets with intermediate busses that have indegree == outdegree == 1 and is not a load bus into a single line\n\nSee remove_bus! for how the two lines are combined.\n\n\n\n\n\n","category":"function"},{"location":"methods/#CommonOPF.trim_tree!","page":"Methods","title":"CommonOPF.trim_tree!","text":"trim_tree!(net::Network)\n\nTrim any branches that have empty busses, i.e. remove the branches that have no loads or DER.\n\n\n\n\n\n","category":"function"},{"location":"methods/#CommonOPF.make_graph","page":"Methods","title":"CommonOPF.make_graph","text":"make_graph(edges::AbstractVector{<:AbstractEdge};  directed::Union{Bool,Missing}=missing)\n\nreturn MetaGraph made up of the edges\n\nAlso the graph[:intbusmap] is created with the dicts for bus => int and int => bus (because Graphs.jl only works with integer nodes)\n\njulia> g[\"13\", :bus]\n10\n\njulia> g[13, :bus]\n\"24\"\n\njulia> get_prop(g, :int_bus_map)[13]\n\"24\"\n\n\n\n\n\n","category":"function"},{"location":"methods/#CommonOPF.leaf_busses","page":"Methods","title":"CommonOPF.leaf_busses","text":"leaf_busses(net::Network)\n\nreturns Vector{String} containing all of the leaf busses in net.graph\n\n\n\n\n\n","category":"function"},{"location":"methods/#BranchFlowModel.solve_metagraph!","page":"Methods","title":"BranchFlowModel.solve_metagraph!","text":"solve_metagraph!(mg::MetaGraphsNext.MetaGraph, builder::Function, tol::T; α::T=0.5, verbose=false) where T <: Real\n\nGiven a MetaGraphsNext.MetaGraph and a JuMP Model builder method iteratively solve the models until the tol is met for the differences provided by BranchFlowModel.get_diffs. \n\nThe builder must accept only one argument of type CommonOPF.AbstractNetwork that returns  a JuMP.AbstractModel. Each model returned from the builder is stored as an :m property in  each vertex of mg.\n\nnote: Note\ntol is compared to the maximum absolute value of all the p, q, and v differences.\n\n\n\n\n\nsolve_metagraph!(mg::MetaGraphsNext.MetaGraph, builder::Function, tols::Vector{T}; α::T=0.5, verbose=false) where T <: Real\n\nGiven a MetaGraphsNext.MetaGraph and a JuMP Model builder method iteratively solve the models until the tols are met for the differences provided by BranchFlowModel.get_diffs. \n\nThe builder must accept only one argument of type CommonOPF.AbstractNetwork that returns  a JuMP.AbstractModel. Each model returned from the builder is stored as an :m property in  each vertex of mg.\n\nnote: Note\nThe tols should have a length of three. The first value is compared to the maximum absolute difference in real power, the second for reactive power, and the third for |v|. All differences are calculated at the leaf/substation connections.\n\n\n\n\n\nsolve_metagraph!(mg::MetaGraphsNext.MetaGraph, builder::Dict{Int64, Function}, tols::Vector{T}; α::T=0.5, verbose=false) where T <: Real\n\nGiven a MetaGraphsNext.MetaGraph and a JuMP Model builder method iteratively solve the models until the tols are met for the differences provided by BranchFlowModel.get_diffs. The builder dict is used to build each model for the corresponding vertex key.\n\nEach function in the builder dict must accept only one argument of type CommonOPF.AbstractNetwork that returns a JuMP.AbstractModel. Each model returned from the builder function is stored as an :m property in each vertex of mg.\n\nnote: Note\nThe tols should have a length of three. The first value is compared to the maximum absolute difference in real power, the second for reactive power, and the third for |v|. All differences are calculated at the leaf/substation connections.\n\n\n\n\n\n","category":"function"},{"location":"methods/#CommonOPF.remove_bus!","page":"Methods","title":"CommonOPF.remove_bus!","text":"remove_bus!(j::String, net::Network{SinglePhase})\n\nRemove bus j in the line i->j->k from the model by making an equivalent line from busses i->k\n\n\n\n\n\nremove_bus!(j::String, net::Network{MultiPhase})\n\nRemove bus j in the line i->j->k from the model by making an equivalent line from busses i->k. We assume the conductors from i->j and j->k have impedance matrices.\n\n\n\n\n\n","category":"function"},{"location":"math/#Single-Phase-BranchFlowModel","page":"Math","title":"Single Phase BranchFlowModel","text":"","category":"section"},{"location":"math/","page":"Math","title":"Math","text":"From [1]","category":"page"},{"location":"math/","page":"Math","title":"Math","text":"Notation:","category":"page"},{"location":"math/","page":"Math","title":"Math","text":"P_ij real power flow from node i to node j\np_j real power injection on node j\n`\\\\mathcal{N}^+ set of all nodes in network except the source\nw_j voltage magnitude squared on node j\nell_ij current magnitude squared  from node i to node j","category":"page"},{"location":"math/","page":"Math","title":"Math","text":"beginaligned\nP_ij - r_ij ell_ij + p_j = sum_kjrightarrow k P_jk  forall j in mathcalN^+ \nQ_ij - x_ij ell_ij + q_j = sum_kjrightarrow k Q_jk  forall j in mathcalN^+ \nw_j = w_i - 2 r_ij P_ij - 2 x_ij Q_ij + (r_ij^2 + x_ij^2) ell_ij  forall j in mathcalN^+ \nw_i ell_ij = P_ij^2 + Q_ij^2 forall (ij) in mathcalE \n(v_jmin)^2 le w_j le (v_jmax)^2  forall j in mathcalN^+ \nendaligned","category":"page"},{"location":"math/#Three-Phase-BranchFlowModel","page":"Math","title":"Three Phase BranchFlowModel","text":"","category":"section"},{"location":"math/","page":"Math","title":"Math","text":"constrain_KVL(m, net::Network{MultiPhase})","category":"page"},{"location":"math/#BranchFlowModel.constrain_KVL-Tuple{Any, Network{MultiPhase}}","page":"Math","title":"BranchFlowModel.constrain_KVL","text":"constrain_KVL(m, net::Network{MultiPhase})\n\nAdd the voltage drop definitions between busses.\n\nw_j = w_i - S_ij Z^star - Z S_ij^star + Z L_ij Z^star\n\n\n\n\n\n","category":"method"},{"location":"math/#References","page":"Math","title":"References","text":"","category":"section"},{"location":"math/#[1]","page":"Math","title":"[1]","text":"","category":"section"},{"location":"math/","page":"Math","title":"Math","text":"Baran, Mesut E., and Felix F. Wu. \"Optimal capacitor placement on radial distribution systems.\" IEEE Transactions on power Delivery 4.1 (1989): 725-734. Chicago\t","category":"page"},{"location":"#BranchFlowModel.jl","page":"User Documentation","title":"BranchFlowModel.jl","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"BranchFlowModel builds the branch flow constraints using JuMP.  The intent of this package is to allow users to build mathematical programs that include BranchFlowModel constraints. No objective is added to the JuMP model in this package and so solving any problem defined by the constraints built by BranchFlowModel.jl is a feasibility problem. Dictionaries of constraints are provided so that one can delete and/or modify the base constraints to fit their problem.","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"warning: Warning\nThis package is under development. Contributions are welcome via fork and pull request.","category":"page"},{"location":"#Inputs","page":"User Documentation","title":"Inputs","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Inputs are defined using CommonOPF.Network structs. ","category":"page"},{"location":"#Building-a-Model","page":"User Documentation","title":"Building a Model","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"build_model!","category":"page"},{"location":"#BranchFlowModel.build_model!","page":"User Documentation","title":"BranchFlowModel.build_model!","text":"build_model!(m::JuMP.AbstractModel, net::Network{SinglePhase})\n\nAdd variables and constraints to m using the values in net. Calls the following functions:\n\nadd_variables(m, net)\nconstrain_power_balance(m, net)\nconstrain_substation_voltage(m, net)\nconstrain_KVL(m, net)\nif relaxed\n    constrain_cone(m, net)\nelse\n    constrain_bilinear(m, net)\nend\n\n\n\n\n\nbuild_model!(m::JuMP.AbstractModel, net::Network{MultiPhase})\n\nAdd variables and constraints to m using the values in net. Calls the following functions:\n\nadd_variables(m, net)\nconstrain_power_balance(m, net)\nconstrain_substation_voltage(m, net)\nconstrain_KVL(m, net)\n\n\n\n\n\n","category":"function"},{"location":"#Variables","page":"User Documentation","title":"Variables","text":"","category":"section"},{"location":"#Single-Phase-Model","page":"User Documentation","title":"Single Phase Model","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Let m be the JuMP.Model provided by the user, then the variables can be accessed via:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"m[:vsqrd] voltage magnitude squared, indexed on busses\nm[:p0], m[:q0] net real, reactive power injection at the substation bus\nm[:Pij], m[:Qij] net real, reactive line flow, indexed on edges\nm[:lij] current magnitude squared, indexed on edges","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"After a model has been solved using JuMP.optimize! variable values can be extracted with JuMP.value. For more see Getting started with JuMP.","category":"page"},{"location":"#MultiPhase-Model","page":"User Documentation","title":"MultiPhase Model","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"The definition of the multiphase variables is done in model_multi_phase.jl as follows:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"# voltage squared is Hermitian\nm[:w] = Dict{Int64, S}()\n# current squared is Hermitian\nm[:l] = Dict{Int64, S}()\n# complex line powers (at the sending end)\nm[:Sij] = Dict{Int64, S}()\n# complex net powers injections \nm[:Sj] = Dict{Int64, S}()\n# Hermitian PSD matrices\nm[:H] = Dict{Int64, S}()","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"where the first key is for the time index and the inner Dict:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"S = Dict{String, AbstractVecOrMat}","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"has string keys for either bus names or edge names.","category":"page"},{"location":"#Accessing-and-Modifying-Constraints","page":"User Documentation","title":"Accessing and Modifying Constraints","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Let the JuMP.Model provided by the user be called m.  Some constraints are stored in the model dict as anonymous constraints with symbol keys.","category":"page"},{"location":"#Power-Injections","page":"User Documentation","title":"Power Injections","text":"","category":"section"},{"location":"","page":"User Documentation","title":"User Documentation","text":"BranchFlowModel.jl uses the convention that power injections are positive (and loads are negative). If no load is provided for a given bus (and phase) then the real and reactive power injections at that bus (and phase) are set to zero with an equality constraint.","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"All power injection constraints are stored in m[:injection_equalities]. The constraints are indexed in the following order:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"by bus name (string), as provided in CommonOPF.busses;\nby \"p\" or \"q\" for real and reactive power respectively;\nby phase number (integer); and\nby time (integer).","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"For example, m[:injection_equalities][\"680\"][\"p\"][2][1] contains the constraint reference for the power injection equality constraint for bus \"680\", real power, in time step 2, on phase 1.","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"If one wished to replace any constraint one must first delete the constraint using the delete function. For example:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"delete(m, m[:cons][:injection_equalities][\"680\"][\"p\"][1])","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"Note that the time index was not provided in the delete command in this example, which implies that the equality constraints for all time steps were deleted. One can also delete individual time step constraints by providing the time index.","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"The deleted constraints can then be replaced with a new set of constraints. For example:","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"m[:cons][:injection_equalities][\"680\"][\"p\"][1] = @constraint(m, [t in 1:net.Ntimesteps],\n    m[:Pj][\"680\",1,t] == -1e3 / net.Sbase\n)","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"where net is the CommonOPF.Network struct for the problem of interest. Note that it is not necessary to store the new constraints in the m[:cons][:injection_equalities].","category":"page"},{"location":"","page":"User Documentation","title":"User Documentation","text":"See the JuMP documentation for more on deleting constraints.","category":"page"}]
}
