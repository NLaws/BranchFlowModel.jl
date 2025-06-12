module BranchFlowModel

using CommonOPF
using JuMP
using LinearAlgebra
import Graphs, MetaGraphsNext
import MetaGraphsNext: vertices
using JSON
import SparseArrays: sparse

# use the same version of MathOptInterface as CommonOPF
# TODO JuMP exports MOI?
const MOI = CommonOPF.MOI
"""
    ModelType

An enum with values:
1. `Unrelaxed`
2. `AngleRelaxation`
3. `Semidefinite`
4. `SOC`  
5. `Linear`

!!! note
    JuMP exports `SecondOrderCone` so we abbreviate it as `SOC`.
"""
@enum ModelType begin
    Unrelaxed
    AngleRelaxation
    Semidefinite
    SOC
    Linear
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
    # recover_voltage_current,  # TODO validate this method
    init_inputs!,
    set_inputs!,
    get_diffs,
    solve_metagraph!,
    check_statuses,
    remove_bus!,

    # ModelType
    Unrelaxed,
    AngleRelaxation,
    Semidefinite,
    SOC,  # JuMP also exports SOC
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
