module BranchFlowModel

using CommonOPF
using JuMP
using LinearAlgebra
import Graphs, MetaGraphsNext
import MetaGraphsNext: vertices
using JSON
import SparseArrays: sparse


"""
    ModelType

An enum with values:
1. `Unrelaxed`
2. `Semidefinite`
3. `SecondOrderCone`
4. `Linear`
"""
@enum ModelType begin
    Unrelaxed
    Semidefinite
    SecondOrderCone
    Linear  # TODO mv LinDistFlow into BranchFlowModel
end


export 
    # new CommonOPF methods
    Network,
    edges,
    busses,
    get_variable_values,
    i_to_j,
    j_to_k,
    rij,
    xij,
    split_at_busses,
    splitting_busses,
    Network_IEEE13_SinglePhase,
    Network_Papavasiliou_2018,
    opf_results,
    reduce_tree!,
    # end new methods TODO distinguish all exports by module

    build_bfm!,
    add_bfm_variables,
    add_linear_variables,
    add_variables,
    add_sdp_variables,
    constrain_power_balance,
    constrain_substation_voltage,
    constrain_KVL,
    constrain_bfm_nlp,
    # recover_voltage_current,  # TODO validate this method
    init_inputs!,
    set_inputs!,
    get_diffs,
    solve_metagraph!,
    check_statuses,
    remove_bus!,

    # ModelType
    Unrelaxed,
    Semidefinite,
    SecondOrderCone,  # JuMP also exports SecondOrderCone
    Linear,

    # MetaGraphsNext
    vertices,
    check_rank_one

include("utils.jl")
include("checks.jl")
include("model_single_phase.jl")
include("model_multi_phase.jl")
include("model_multi_phase_linear.jl")
include("decomposition.jl")

end
