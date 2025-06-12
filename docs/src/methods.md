# Methods
Some various methods used in BranchFlowModel.jl:

```@contents
Pages = ["methods.md"]
Depth = 3
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
BranchFlowModel.add_vsqrd_variables
BranchFlowModel.add_isqrd_variables
```

## Constraint builders
```@docs
BranchFlowModel.constrain_power_balance
BranchFlowModel.constrain_KVL
BranchFlowModel.constrain_bfm_nlp
BranchFlowModel.constrain_linear_power_balance
BranchFlowModel.constrain_KVL_linear
BranchFlowModel.constrain_power_balance_with_isqrd_losses
BranchFlowModel.constrain_bilinear
BranchFlowModel.constrain_substation_voltage
BranchFlowModel.constrain_cone
```

## Other
```@docs
BranchFlowModel.MPij
BranchFlowModel.MQij
check_rank_one
reduce_tree!
trim_tree!
set_inputs!
get_diffs
solve_metagraph!
check_statuses
remove_bus!
```