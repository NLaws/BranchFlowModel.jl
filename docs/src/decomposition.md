Methods to decompose and solve the `SinglePhase` Branch Flow Model are provided based on the work in 
[[2]](@ref). These methods are most advantageous when solving the non-linear (unrelaxed) power flow
equations and are only valid(?) in radial networks.

```@docs
init_inputs!
set_inputs!
get_diffs
split_inputs
```

### [2]
Sadnan, Rabayet, and Anamika Dubey. "Distributed optimization using reduced network equivalents for radial power distribution systems." IEEE Transactions on Power Systems 36.4 (2021): 3645-3656.