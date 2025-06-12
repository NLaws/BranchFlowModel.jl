# Single Phase Models
BranchFlowModel.jl provides methods to build many different variations of the Branch Flow Model,
including single phase and multiphase models. Each of the model types supported are documented below.
```@contents
Pages = ["single_phase_models.md"]
Depth = 3
```
```@setup imports
using BranchFlowModel
using CommonOPF
using JuMP
```


## `AngleRelaxation`
The `Unrelaxed` multiphase model is build by passing a `JuMP.Model`, `Network{SinglePhase}`, and the
`AngleRelaxation` type to [`build_bfm!`](@ref).

```@example imports
net = CommonOPF.Network_IEEE13_SinglePhase()
m = JuMP.Model()

build_bfm!(m, net, AngleRelaxation)
println("Variable information:")
CommonOPF.print_var_info(net)
println("Constraint information:")
CommonOPF.print_constraint_info(net)
```

The [`build_bfm!](@ref) method uses:
- [`BranchFlowModel.add_linear_variables`](@ref)
- [`BranchFlowModel.add_vsqrd_variables`](@ref)
- [`BranchFlowModel.add_isqrd_variables`](@ref)
- [`BranchFlowModel.constrain_power_balance_with_isqrd_losses`](@ref)
- [`BranchFlowModel.constrain_KVL`](@ref)
- [`BranchFlowModel.constrain_bilinear`](@ref)

The math underlying the model is as follows [[Farivar and Low]](@ref):

```math
\begin{aligned}
    &\sum_{i : i \rightarrow j} ( p_{ij} - r_{ij} \ell_{ij} ) + p_j 
    = \sum_{k : j \rightarrow k} p_{jk} 
    \quad \forall j \in \mathcal{N}
    \\
    &\sum_{i : i \rightarrow j} ( q_{ij} - x_{ij} \ell_{ij} ) + q_j 
    = \sum_{k : j \rightarrow k} q_{jk} 
    \quad \forall j \in \mathcal{N}
    \\
    &w_{j} = w_{i} - 2 (p_{ij} r_{ij} + q_{ij} x_{ij}) + (r_{ij}^2 + x_{ij}^2) \ell_{ij} 
    \quad \forall (i, j) \in \mathcal{E}
    \\
    & \ell_{ij} w_i = p_{ij}^2 + q_{ij}^2 \quad \forall (i, j) \in \mathcal{E}
\end{aligned}
```


## `SecondOrderCone` models
```math
\begin{aligned}
    &\sum_{i : i \rightarrow j} ( p_{ij} - r_{ij} \ell_{ij} ) + p_j 
    = \sum_{k : j \rightarrow k} p_{jk} 
    \quad \forall j \in \mathcal{N}
    \\
    &\sum_{i : i \rightarrow j} ( q_{ij} - x_{ij} \ell_{ij} ) + q_j 
    = \sum_{k : j \rightarrow k} q_{jk} 
    \quad \forall j \in \mathcal{N}
    \\
    &w_{j} = w_{i} - 2 (p_{ij} r_{ij} + q_{ij} x_{ij}) + (r_{ij}^2 + x_{ij}^2) \ell_{ij} 
    \quad \forall (i, j) \in \mathcal{E}
    \\
    & \ell_{ij} \geq \frac{p_{ij}^2 + q_{ij}^2}{w_i} \quad \forall (i, j) \in \mathcal{E}
\end{aligned}
```


## `Linear` models
The single phase "LinDistFlow" model from [[Baran and Wu]](@ref)

Notation:
- ``P_{ij}`` real power flow from node ``i`` to node ``j``
- ``p_j`` real power injection on node ``j``
- ``\mathcal{N}^+`` set of all nodes in network except the source
- ``w_j`` voltage magnitude squared on node ``j``

```math
\begin{aligned}
P_{ij} + p_j = \sum_{k:j\rightarrow k} P_{jk} \ \forall j \in \mathcal{N}^+ \\
Q_{ij} + q_j = \sum_{k:j\rightarrow k} Q_{jk} \ \forall j \in \mathcal{N}^+ \\
w_j = w_i - 2 r_{ij} P_{ij} - 2 x_{ij} Q_{ij} \ \forall j \in \mathcal{N}^+ \\
(v_{j,\min})^2 \le w_j \le (v_{j,\max})^2 \ \forall j \in \mathcal{N}^+ 
\end{aligned}
```

# References

### [Baran and Wu]
Baran, Mesut E., and Felix F. Wu. "Optimal capacitor placement on radial distribution systems." IEEE Transactions on power Delivery 4.1 (1989): 725-734.
Chicago	

### [Farivar and Low]
Farivar, Masoud, and Steven H. Low. "Branch flow model: Relaxations and convexification—Part I." IEEE Transactions on Power Systems 28.3 (2013): 2554-2564.
