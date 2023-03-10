module BranchFlowModel

using JuMP
using LinearAlgebra
import SparseArrays: sparse
import PowerModelsDistribution: parse_dss, DELTA

const DEFAULT_AMP_LIMIT = 1000.0

export 
    Inputs,
    singlephase38linesInputs,
    dsstxt_to_sparse_array, 
    parse_dss,
    build_model!,
    add_variables,
    constrain_power_balance,
    constrain_substation_voltage,
    constrain_KVL,
    constrain_loads,
    constrain_bounds,
    get_bus_values, 
    get_edge_values, 
    get_ijlinecode,
    # recover_voltage_current,  # TODO validate this method
    i_to_j, 
    j_to_k, 
    rij, 
    xij,
    zij,
    check_soc_inequalities,
    get_load_bal_shadow_prices

include("types.jl")
include("io.jl")
include("inputs.jl")
include("utils.jl")
include("model_single_phase.jl")
include("model_multi_phase.jl")

end
