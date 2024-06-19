
# Single Phase BranchFlowModel
From [[1]](@ref)

Notation:
- ``P_{ij}`` real power flow from node ``i`` to node ``j``
- ``p_j`` real power injection on node ``j``
- ``\\mathcal{N}^+` set of all nodes in network except the source
- ``w_j`` voltage magnitude squared on node ``j``
- ``\ell_{ij}`` current magnitude squared  from node ``i`` to node ``j``

```math
\begin{aligned}
P_{ij} - r_{ij} \ell_{ij} + p_j = \sum_{k:j\rightarrow k} P_{jk} \ \forall j \in \mathcal{N}^+ \\
Q_{ij} - x_{ij} \ell_{ij} + q_j = \sum_{k:j\rightarrow k} Q_{jk} \ \forall j \in \mathcal{N}^+ \\
w_j = w_i - 2 r_{ij} P_{ij} - 2 x_{ij} Q_{ij} + (r_{ij}^2 + x_{ij}^2) \ell_{ij} \ \forall j \in \mathcal{N}^+ \\
w_i \ell_{ij} = P_{ij}^2 + Q_{ij}^2 \forall (i,j) \in \mathcal{E} \\
(v_{j,\min})^2 \le w_j \le (v_{j,\max})^2 \ \forall j \in \mathcal{N}^+ 
\end{aligned}
```

# Three Phase BranchFlowModel
```@docs
constrain_KVL(m, net::Network{MultiPhase})
```


# References

### [1]
Baran, Mesut E., and Felix F. Wu. "Optimal capacitor placement on radial distribution systems." IEEE Transactions on power Delivery 4.1 (1989): 725-734.
Chicago	
