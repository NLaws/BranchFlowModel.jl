module BranchFlowModel

using CommonOPF
import CommonOPF: rij, xij  # extending to MultiPhase, need to consolidate into CommonOPF
using JuMP
using LinearAlgebra
using Graphs, MetaGraphs
using JSON
import SparseArrays: sparse

const DEFAULT_AMP_LIMIT = 1000.0

export 
    Inputs,
    # new CommonOPF methods
    Network,
    edges,
    busses,
    i_to_j,
    j_to_k,
    rij, 
    xij,
    Network_IEEE13_SinglePhase,
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
    check_rank_one,
    get_variable_values,
    get_edge_values, 
    get_ijlinecode,
    get_ijlinelength,
    get_ij_idx,
    # recover_voltage_current,  # TODO validate this method
    zij,
    check_soc_inequalities,
    check_connected_graph,
    get_load_bal_shadow_prices,
    voltage_values_by_time_bus,
    current_values_by_time_edge,
    line_flow_values_by_time_edge,
    reduce_tree!,
    make_graph,
    init_inputs!,
    set_inputs!,
    split_at_busses,
    get_diffs,
    all_outneighbors,
    all_inneighbors,
    busses_from_deepest_to_source,
    vertices_from_deepest_to_source,
    splitting_busses,
    connect_subgraphs_at_busses,
    split_inputs,
    build_metagraph,
    solve_metagraph!,
    metagraph_voltages,
    check_unique_solution_conditions,
    check_statuses,
    reg_busses,
    turn_ratio,
    vreg,
    has_vreg,
    combine_parallel_lines!,
    remove_bus!,
    paths_between,
    Results

include("inputs.jl")
include("utils.jl")
include("results.jl")
include("checks.jl")
include("model_single_phase.jl")
include("model_single_phase_network.jl")
include("model_multi_phase.jl")
include("decomposition.jl")

end
