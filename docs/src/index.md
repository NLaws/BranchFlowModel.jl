# BranchFlowModel.jl

`BranchFlowModel` builds the branch flow constraints using JuMP. The intent of this package is to
allow users to build mathematical programs that include BranchFlowModel constraints. No objective is
added to the JuMP model in this package and so solving any problem defined by the constraints built
by BranchFlowModel.jl is a feasibility problem. Dictionaries of constraints are provided so that one
can delete and/or modify the base constraints to fit their problem.


!!! warning
    This package is under development. Contributions are welcome via fork and pull request.

# Inputs
Inputs are defined using [`CommonOPF.Network` structs](https://nlaws.github.io/CommonOPF.jl/stable/network/). 


# Building a Model
Building a `BranchFlowModel` requires three things:
1. a JuMP Model,
2. a `CommonOPF.Network`, and
3. the type of model to be built, i.e. one of the [`BranchFlowModel.ModelType`](@ref)
```@docs
BranchFlowModel.ModelType
```
To build a model see [build_bfm!](@ref)

