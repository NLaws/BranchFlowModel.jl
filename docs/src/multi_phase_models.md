# Multiphase Models
BranchFlowModel.jl provides methods to build many different variations of the Branch Flow Model,
including single phase and multiphase models. Each of the model types supported are documented below.
```@contents
Pages = ["multi_phase_models.md"]
Depth = 3
```
```@setup imports
using BranchFlowModel
using CommonOPF
using JuMP
```


## `Unrelaxed` models
The `Unrelaxed` multiphase model is build by passing a `JuMP.Model`, `Network{MultiPhase}`, and the
`Unrelaxed` type to [`build_bfm!`](@ref).

```@example imports
net = CommonOPF.Network_IEEE13()
m = JuMP.Model()

build_bfm!(m, net, Unrelaxed)
println("Variable information:")
CommonOPF.print_var_info(net)
println("Constraint information:")
CommonOPF.print_constraint_info(net)
```

The [`build_bfm!](@ref) method uses:
- [`add_bfm_variables`](@ref)
- [`constrain_bfm_nlp`](@ref)
- [`constrain_power_balance`](@ref)

The math underlying the model is as follows:
```math
\begin{aligned}
    &\boldsymbol S_{ij} = \boldsymbol v_i^{\Phi_{ij}} \boldsymbol i_{ij}^H
    \quad \forall (i, j) \in \mathcal{E}
    \\
    &\boldsymbol v_i^{\Phi_{ij}} - \boldsymbol v_j = \boldsymbol Z_{ij} \boldsymbol i_{ij}
    \quad \forall (i, j) \in \mathcal{E}
    \\
    &\sum_{i : i \rightarrow j}  \text{diag}( \boldsymbol S_{ij} - \boldsymbol Z_{ij} \left[ \boldsymbol i_{ij} \boldsymbol i_{ij}^H \right]) 
    + \boldsymbol s_j 
    = \sum_{k : j \rightarrow k} \text{diag}( \boldsymbol S_{jk} )^{\Phi_j}
    \quad \forall j \in \mathcal{N}
\end{aligned}
```
For the nomenclatures see TODO.


## `Semidefinite` models
The `Semidefinite` multiphase model is build by passing a `JuMP.Model`, `Network{MultiPhase}`, and the
`Semidefinite` type to [`build_bfm!`](@ref).


```@example imports
net = CommonOPF.Network_IEEE13()
m = JuMP.Model()

build_bfm!(m, net, Semidefinite)
println("Variable information:")
CommonOPF.print_var_info(net)
println("Constraint information:")
CommonOPF.print_constraint_info(net)
```

The [`build_bfm!](@ref) method uses:
- [`add_sdp_variables`](@ref)
- [`constrain_KVL`](@ref)
- [`constrain_power_balance`](@ref)

## `Linear` models
TODO document constraints and their accessors as well. 
```@example imports
net = CommonOPF.Network_IEEE13()
m = JuMP.Model()

build_bfm!(m, net, Linear)
CommonOPF.print_var_info(net)
```
