# Methods
Some various methods used in BranchFlowModel.jl:

!!! warning
    This list of exported methods may not be up to date and there are missing doc strings.
    Contributions are welcome via fork and pull request.

```@docs
Inputs(::String, ::String)
Inputs(::AbstractVector{<:Tuple}, ::AbstractVector{<:AbstractString}, ::AbstractVector{<:Real}, ::AbstractVector{<:AbstractVector}, ::String)
singlephase38linesInputs
dsstxt_to_sparse_array 
dss_files_to_dict
build_model!(m::JuMP.AbstractModel, p::Inputs{BranchFlowModel.SinglePhase})
build_model!(m::JuMP.AbstractModel, p::Inputs{BranchFlowModel.MultiPhase})
add_variables
constrain_power_balance
constrain_substation_voltage
constrain_KVL
constrain_loads
constrain_bounds
check_rank_one
get_bus_values 
get_edge_values 
get_ijlinecode
get_ijlinelength
get_ij_idx
i_to_j 
j_to_k 
rij 
xij
zij
check_soc_inequalities
check_connected_graph
get_load_bal_shadow_prices
voltage_values_by_time_bus
current_values_by_time_edge
line_flow_values_by_time_edge
reduce_tree!
trim_tree!
make_graph
leaf_busses
init_inputs!
set_inputs!
split_at_busses
get_diffs
splitting_busses
connect_subgraphs_at_busses
split_inputs
build_metagraph
solve_metagraph!
metagraph_voltages
check_unique_solution_conditions
check_statuses
reg_busses
turn_ratio
vreg
has_vreg
remove_bus!
paths_between
Results(m::JuMP.AbstractModel, p::Inputs{SinglePhase}; digits=8)
```