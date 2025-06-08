# Variables
BranchFlowModel.jl documents each model's variables using the `CommonOPF.VariableInfo`, stored in the
`CommonOPF.Network.var_info` dict. The `keys(var_info)` are symbols that correspond with the symbols
stored in the `JuMP.Model.obj_dict`. If `m` is the `JuMP.Model` provided by the user, then the
variables can be accessed via `m[:a_var_symbol]`. Similarly, the `VariableInfo` is accessed via `network.var_info[:a_var_symbol]`.

After a model has been solved using `JuMP.optimize!` variable values can be extracted with
`JuMP.value`. For more see [Getting started with
JuMP](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/getting_started_with_JuMP/#Getting-started-with-JuMP).

Variables are stored in dictionaries. The order of the keys in the dictionaries is standardized in
`CommonOPF` variable containers.
