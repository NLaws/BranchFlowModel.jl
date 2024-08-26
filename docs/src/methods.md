# Methods
Some various methods used in BranchFlowModel.jl:

!!! warning
    This list of exported methods may not be up to date and there are missing doc strings.
    Contributions are welcome via fork and pull request.

```@docs
build_bfm!
constrain_power_balance
constrain_KVL
check_rank_one
check_soc_inequalities
get_load_bal_shadow_prices
reduce_tree!
trim_tree!
set_inputs!
get_diffs
solve_metagraph!
metagraph_voltages
check_statuses
remove_bus!
```