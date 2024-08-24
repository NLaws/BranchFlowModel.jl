module BranchFlowModel

using CommonOPF
using JuMP
using LinearAlgebra
import Graphs, MetaGraphsNext
import MetaGraphsNext: vertices
using JSON
import SparseArrays: sparse


@enum ModelType begin
    Unrelaxed
    Semidefinite
    SecondOrderCone
    # Linear  # TODO mv LinDistFlow into BranchFlowModel
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
    split_at_busses,
    splitting_busses,
    xij,
    Network_IEEE13_SinglePhase,
    Network_Papavasiliou_2018,
    Results,
    reduce_tree!,
    # end new methods TODO distinguish all exports by module

    build_bfm!,
    add_variables,
    constrain_power_balance,
    constrain_substation_voltage,
    constrain_KVL,
    constrain_bounds,
    # recover_voltage_current,  # TODO validate this method
    make_graph,
    init_inputs!,
    set_inputs!,
    get_diffs,
    solve_metagraph!,
    metagraph_voltages,
    check_statuses,
    reg_busses,
    remove_bus!,

    # ModelType
    Unrelaxed,
    Semidefinite,
    SecondOrderCone,  # JuMP also exports SecondOrderCone

    # MetaGraphsNext
    vertices,
    check_rank_one

include("inputs.jl")
include("utils.jl")
include("checks.jl")
include("model_single_phase.jl")
include("model_multi_phase.jl")
include("decomposition.jl")

end
