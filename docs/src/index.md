# BranchFlowModel.jl

`BranchFlowModel` builds the branch flow constraints using JuMP. 
The intent of this package is to allow users to build mathematical programs that include BranchFlowModel constraints.
No objective is added to the JuMP model in this package and so solving any problem defined by the constraints built by BranchFlowModel.jl is a feasibility problem.
Dictionaries of constraints are provided so that one can delete and/or modify the base constraints to fit their problem.


!!! warning
    This package is under development. Contributions are welcome via fork and pull request.

# Inputs
There are two methods for creating `Inputs`:
1. Using openDSS files
2. Providing the network topology
```@docs
Inputs(::String, ::String)
Inputs(::AbstractVector{<:Tuple}, ::AbstractVector{<:AbstractString}, ::AbstractVector{<:Real}, ::AbstractVector{<:AbstractVector}, ::String)
```
Both of the `Inputs` functions return a mutable `Inputs` struct:
```@docs
Inputs
```

# Building a Model
The `build_ldf!` function takes a `JuMP.Model` and `Inputs` struct as its two arguments and adds the variables and constraints:
```@docs
build_model!
```

# Variables

## Single Phase Model
Let `m` be the JuMP.Model provided by the user, then the variables can be accessed via:
- `m[:vsqrd]` voltage magnitude squared, indexed on busses, time
- `m[:Pj], m[:Qj]` net real, reactive power injection, indexed on busses, time
- `m[:Pij], m[:Qij]` net real, reactive line flow, indexed on edges, time
- `m[:lij]` current magnitude squared, indexed on edges, time
Note that the `m[:Pj], m[:Qj]` are not truly variables since they are defined by the loads (unless you modify the model to make them decisions).
After a model has been solved using `JuMP.optimize!` variable values can be extracted with `JuMP.value`. For more see [Getting started with JuMP](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/getting_started_with_JuMP/#Getting-started-with-JuMP).

## MultiPhase Model
The definition of the multiphase variables is done in `model_multi_phase.jl` as follows:
```julia
# voltage squared is Hermitian
m[:w] = Dict{Int64, S}()
# current squared is Hermitian
m[:l] = Dict{Int64, S}()
# complex line powers (at the sending end)
m[:Sij] = Dict{Int64, S}()
# complex net powers injections 
m[:Sj] = Dict{Int64, S}()
# Hermitian PSD matrices
m[:H] = Dict{Int64, S}()
```
where the first key is for the time index and the inner `Dict`:
```julia
S = Dict{String, AbstractVecOrMat}
```
has string keys for either bus names or edge names, (which are stored in the `Inputs` as `Inputs.busses` and `Inputs.edge_keys` respectively).

# Accessing and Modifying Constraints
Let the JuMP.Model provided by the user be called `m`. 
Some constraints are stored in the model dict as anonymous constraints with symbol keys.

## Power Injections
BranchFlowModel.jl uses the convention that power injections are positive (and loads are negative). If no load is provided for a given bus (and phase) then the real and reactive power injections at that bus (and phase) are set to zero with an equality constraint.

All power injection constraints are stored in `m[:injection_equalities]`. The constraints are indexed in the following order:
1. by bus name (string), as provided in `Inputs.busses`;
2. by `"p"` or `"q"` for real and reactive power respectively;
3. by phase number (integer); and
4. by time (integer).
For example, `m[:injection_equalities]["680"]["p"][2][1]` contains the constraint reference for the power injection equality constraint for bus "680", real power, in time step 2, on phase 1.

If one wished to replace any constraint one must first delete the constraint using the `delete` function. For example:
```julia
delete(m, m[:cons][:injection_equalities]["680"]["p"][1])
```
Note that the time index was not provided in the `delete` command in this example, which implies that the equality constraints for all time steps were deleted. One can also delete individual time step constraints by providing the time index.

The deleted constraints can then be replaced with a new set of constraints. For example:
```julia
m[:cons][:injection_equalities]["680"]["p"][1] = @constraint(m, [t in 1:p.Ntimesteps],
    m[:Pj]["680",1,t] == -1e3 / p.Sbase
)
```
where `p` is short for "parameters" and is the `Inputs` struct for the problem of interest. Note that it is not necessary to store the new constraints in the `m[:cons][:injection_equalities]`.

See the [JuMP documentation](https://jump.dev/JuMP.jl/stable/manual/constraints/#Delete-a-constraint) for more on deleting constraints.

# Results
```@docs
Results(m::AbstractModel, p::Inputs{SinglePhase})
```
