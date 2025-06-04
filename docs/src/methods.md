# Methods
Some various methods used in BranchFlowModel.jl:

!!! warning
    This list of exported methods may not be up to date and there are missing doc strings.
    Contributions are welcome via fork and pull request.

## Model builders
```@docs
build_bfm!
```

## Variable builders
```@docs
add_bfm_variables
add_sdp_variables
add_linear_variables
```

## Constraint builders
```@docs
constrain_power_balance
constrain_KVL
```

## Other
```@docs
check_rank_one
reduce_tree!
trim_tree!
set_inputs!
get_diffs
solve_metagraph!
check_statuses
remove_bus!
```