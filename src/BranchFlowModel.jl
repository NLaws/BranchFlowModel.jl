module BranchFlowModel

using CommonOPF
import CommonOPF: rij, xij  # extending to MultiPhase, need to consolidate into CommonOPF
using JuMP
using LinearAlgebra
import Graphs, MetaGraphsNext
using JSON
import SparseArrays: sparse

const DEFAULT_AMP_LIMIT = 1000.0

export 
    Inputs,
    # new CommonOPF methods
    Network,
    edges,
    busses,
    get_variable_values,
    i_to_j,
    j_to_k,
    rij, 
    xij,
    Network_IEEE13_SinglePhase,
    Network_Papavasiliou_2018,
    # end new methods TODO distinguish all exports by module
    singlephase38linesInputs,
    dsstxt_to_sparse_array, 
    dss_files_to_dict,
    build_model!,
    add_variables,
    constrain_power_balance,
    constrain_substation_voltage,
    constrain_KVL,
    constrain_loads,
    constrain_bounds,
    get_edge_values, 
    get_ijlinecode,
    get_ijlinelength,
    get_ij_idx,
    # recover_voltage_current,  # TODO validate this method
    zij,
    voltage_values_by_time_bus,
    current_values_by_time_edge,
    line_flow_values_by_time_edge,
    reduce_tree!,
    make_graph,
    init_inputs!,
    set_inputs!,
    get_diffs,
    split_inputs,
    solve_metagraph!,
    metagraph_voltages,
    check_statuses,
    reg_busses,
    turn_ratio,
    vreg,
    has_vreg,
    remove_bus!,
    paths_between,
    Results

include("inputs.jl")
include("utils.jl")
include("checks.jl")
include("model_single_phase_network.jl")
include("decomposition.jl")

end
