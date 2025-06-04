# Multiphase Models
BranchFlowModel.jl provides methods to build many different variations of the Branch Flow Model,
including single phase and multiphase models. Each of the model types supported are documented below.
```@contents
Pages = ["multi_phase_models.md"]
Depth = 3
```

The definition of the multiphase variables is done in `model_multi_phase.jl` as follows:

```@setup print_var_info-multiphase
# make the most basic multiphase model so that we can add variables and print their info
using BranchFlowModel
using CommonOPF
using JuMP


net = Network(Dict(
    :Network => Dict(
        :substation_bus => "source"
    ),
    :Conductor => [Dict(
        :busses => ("source", "b1"),
        :phases => [1,2,3]
    )]
))

```


## `Unrelaxed` models
```@docs
add_bfm_variables
```


## `Semidefinite` models
```@docs
add_sdp_variables
```

## `Linear` models
TODO use CommonOPF test networks and `build_bfm` instead of `add_linear_variables` to show how to
build models too. Eventually document constraints and their accessors as well. Maybe organize docs
by model type rather than variables and constraints -- or make the variables and constraints
sections general (with some specific examples), and have the model type documentation with repeated
format of docs.
```@example print_var_info-multiphase
m = JuMP.Model()
println(typeof(net))
add_linear_variables(m, net)
CommonOPF.print_var_info(net)
```
