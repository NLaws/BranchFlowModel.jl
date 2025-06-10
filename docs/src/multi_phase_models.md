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
- [`BranchFlowModel.add_bfm_variables`](@ref)
- [`BranchFlowModel.constrain_bfm_nlp`](@ref)
- [`BranchFlowModel.constrain_power_balance`](@ref)

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
For the nomenclature see TODO.


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
- [`BranchFlowModel.add_sdp_variables`](@ref)
- [`BranchFlowModel.constrain_KVL`](@ref)
- [`BranchFlowModel.constrain_power_balance`](@ref)

## `Linear` models
The `Linear` multiphase model is build by passing a `JuMP.Model`, `Network{MultiPhase}`, and the
`Linear` type to [`build_bfm!`](@ref).

```@example imports
net = CommonOPF.Network_IEEE13()
m = JuMP.Model()

build_bfm!(m, net, Linear)
println("Variable information:")
CommonOPF.print_var_info(net)
println("Constraint information:")
CommonOPF.print_constraint_info(net)
```

The [`build_bfm!](@ref) method uses:
- [`BranchFlowModel.add_linear_variables`](@ref)
- [`BranchFlowModel.constrain_linear_power_balance`](@ref)
- [`BranchFlowModel.constrain_KVL_linear`](@ref)


The math underlying the model is as follows (from [Arnold et al.](@ref)):
```math
\begin{aligned}
P_{ij,\phi} + p_{j,\phi} = \sum_{k:j\rightarrow k} P_{jk,\phi} \ \forall j \in \mathcal{N}^+, \forall \phi \in [1,2,3] \\
Q_{ij,\phi} + q_{j,\phi} = \sum_{k:j\rightarrow k} Q_{jk,\phi} \ \forall j \in \mathcal{N}^+, \forall \phi \in [1,2,3] \\
\boldsymbol{w}_j = \boldsymbol{w}_i + \boldsymbol{M}_{P,ij} \boldsymbol{P}_{ij} + \boldsymbol{M}_{Q,ij} \boldsymbol{Q}_{ij} \\
(\boldsymbol{v}_{j,\min})^2 \le \boldsymbol{w}_j \le (\boldsymbol{v}_{j,\max})^2 \ \forall j \in \mathcal{N}^+ \\
\boldsymbol{M}_{P,ij} = \begin{bmatrix}
-2r_{11}                & r_{12}-\sqrt{3}x_{12} & r_{13}+\sqrt{3}x_{13} \\
  r_{21}+\sqrt{3}x_{21} & -2r_{22} & r_{23}-\sqrt{3}x_{23} \\
  r_{31}-\sqrt{3}x_{31} & r_{32}+\sqrt{3}x_{32} & -2r_{33}
\end{bmatrix} \\
\boldsymbol{M}_{Q,ij} = \begin{bmatrix}
-2x_{11}                &   x_{12}+\sqrt{3}r_{12} &   x_{13}-\sqrt{3}r_{13} \\
  x_{21}-\sqrt{3}r_{21}  & -2x_{22}                &   x_{23}+\sqrt{3}r_{23} \\
  x_{31}+\sqrt{3}r_{31} &   x_{32}-\sqrt{3}r_{32} & -2x_{33}
\end{bmatrix} 
\end{aligned}
```

# References

### Arnold et al.
Arnold, Daniel B., et al. "Optimal dispatch of reactive power for voltage regulation and balancing in unbalanced distribution systems." 2016 IEEE Power and Energy Society General Meeting (PESGM). IEEE, 2016.