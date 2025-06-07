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

```@example imports
net = CommonOPF.Network_IEEE13()
m = JuMP.Model()

build_bfm!(m, net, Unrelaxed)
CommonOPF.print_var_info(net)
```


## `Semidefinite` models

```@example imports
net = CommonOPF.Network_IEEE13()
m = JuMP.Model()

build_bfm!(m, net, Semidefinite)
CommonOPF.print_var_info(net)
```

## `Linear` models
TODO document constraints and their accessors as well. 
```@example imports
net = CommonOPF.Network_IEEE13()
m = JuMP.Model()

build_bfm!(m, net, Linear)
CommonOPF.print_var_info(net)
```
