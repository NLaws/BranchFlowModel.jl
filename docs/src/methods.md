# Methods
Some various methods used in BranchFlowModel.jl:

```@contents
Pages = ["methods.md"]
Depth = 2
```

## Model builders
```@docs
build_bfm!
```

## Variable builders
```@docs
BranchFlowModel.add_bfm_variables
BranchFlowModel.add_sdp_variables
BranchFlowModel.add_linear_variables
```

## Constraint builders
```@docs
BranchFlowModel.constrain_power_balance
BranchFlowModel.constrain_KVL
BranchFlowModel.constrain_bfm_nlp
BranchFlowModel.constrain_linear_power_balance
BranchFlowModel.constrain_KVL_linear
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