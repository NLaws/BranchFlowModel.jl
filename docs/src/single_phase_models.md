# Single Phase Models
BranchFlowModel.jl provides methods to build many different variations of the Branch Flow Model,
including single phase and multiphase models. Each of the model types supported are documented below.
```@contents
Pages = ["single_phase_models.md"]
Depth = 3
```

## Angle relaxation
Let `m` be the JuMP.Model provided by the user, then the variables can be accessed via:
- `m[:vsqrd]` voltage magnitude squared, indexed on busses
- `m[:p0], m[:q0]` net real, reactive power injection at the substation bus
- `m[:pij], m[:qij]` net real, reactive line flow, indexed on edges
- `m[:lij]` current magnitude squared, indexed on edges


## `Linear` models
