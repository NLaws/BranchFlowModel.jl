# Methods
Some various methods used in BranchFlowModel.jl:

!!! warning
    This list of exported methods may not be up to date and there are missing doc strings.
    Contributions are welcome via fork and pull request.

```@docs
dsstxt_to_sparse_array 
dss_files_to_dict
build_model!(m::JuMP.AbstractModel, net::Network{SinglePhase}; relaxed::Bool=true)build_model!(m::JuMP.AbstractModel, net::Network{MultiPhase}; PSD::Bool=true)add_variables
constrain_power_balance
constrain_substation_voltage
constrain_KVL
constrain_bounds
check_rank_one
get_bus_values 
get_edge_values 
check_soc_inequalities
get_load_bal_shadow_prices
current_values_by_time_edge
line_flow_values_by_time_edge
reduce_tree!
trim_tree!
make_graph
leaf_busses
set_inputs!
get_diffs
solve_metagraph!
metagraph_voltages
check_unique_solution_conditions
check_statuses
reg_busses
remove_bus!
```