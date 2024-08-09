# Methods
Some various methods used in BranchFlowModel.jl:

!!! warning
    This list of exported methods may not be up to date and there are missing doc strings.
    Contributions are welcome via fork and pull request.

```@docs
build_model!
constrain_power_balance
constrain_substation_voltage
constrain_KVL
constrain_bounds
check_rank_one
get_bus_values 
check_soc_inequalities
get_load_bal_shadow_prices
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