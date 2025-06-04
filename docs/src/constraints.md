
# Accessing and Modifying Constraints
Let the JuMP.Model provided by the user be called `m`. 
Some constraints are stored in the model dict as anonymous constraints with symbol keys.

## Power Injections
BranchFlowModel.jl uses the convention that power injections are positive (and loads are negative). If no load is provided for a given bus (and phase) then the real and reactive power injections at that bus (and phase) are set to zero with an equality constraint.

All power injection constraints are stored in `m[:injection_equalities]`. The constraints are indexed in the following order:
1. by bus name (string), as provided in `CommonOPF.busses`;
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
m[:cons][:injection_equalities]["680"]["p"][1] = @constraint(m, [t in 1:net.Ntimesteps],
    m[:Pj]["680",1,t] == -1e3 / net.Sbase
)
```
where `net` is the `CommonOPF.Network` struct for the problem of interest. Note that it is not necessary to store the new constraints in the `m[:cons][:injection_equalities]`.

See the [JuMP documentation](https://jump.dev/JuMP.jl/stable/manual/constraints/#Delete-a-constraint) for more on deleting constraints.
