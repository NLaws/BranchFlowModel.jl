# BranchFlowModel Changelog

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
