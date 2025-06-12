# BranchFlowModel Changelog

## dev
- upgrade CommonOPF to v0.4.6
- refactor `SecondOrderCone` -> `SOC` to avoid collision with `JuMP.SecondOrderCone`
- SinglePhase `Unrelaxed` model is now `AngleRelaxation` (as it should be)
- ported the multiphase LinDistFlow validation test from LinDistFlow.jl
- refactored and added symbols used to store constraint references in the `JuMP.Model`
    - `:kvl` -> `:kvl_constraints`
        - these were erroneously indexed on busses instead of edges
    - `:loadbalcons` -> `:power_balance_constraints`
        - also changed the `"p"` and `"q"` keys to `:real` and `:reactive`
    - add `:sending_power_flow_constraints`
- added `net.var_info` and `net.constraint_info` for multiphase models and improved docs
- stopped exporting the variable and constraint builder methods

## v0.4.5
- move some model building utilities to CommonOPF for reuse in BusInjectionModel

## v0.4.4
- transition to CommonOPF.Network (from CommonOPF.Inputs) in v0.4 of CommonOPF

## v0.4.3
- upgrade CommonOPF to v0.3.8

## v0.4.2
- upgrade CommonOPF to v0.3.7

## v0.4.1
- mv `remove_bus!` and `reduce_tree!` to CommonOPF
- mv SinglePhase `rij` and `xij` to CommonOPF

## v0.4.0
- deprecate `get_bus_values`
- mv `get_variable_values` to CommonOPF

## v0.3.2
- add `shunt_susceptance` times `vsqrd` to flow balance equations in single phase model
- add `Results.shadow_prices`

## v0.3.1
- change JuMP compat to "1" (was "1.9") for better cross-compatibility

## v0.3.0
- CommonOPF upgrade to 0.3

## v0.2.1
- add `Results`
- improve `combine_parallel_lines!`
- mv `reg_busses`, `turn_ratio`, `has_vreg`, `vreg` to CommonOPF

## v0.2.0
- move io.jl, types.jl, inputs.jl and some utils.jl to CommonOPF (new dependency)

## v0.1.0
- initial release
