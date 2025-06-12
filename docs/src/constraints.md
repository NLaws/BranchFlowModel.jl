
# Accessing and Modifying Constraints
Let the JuMP.Model provided by the user be called `m`. 
The constraints are stored in the model dict as anonymous constraints with symbol keys.
The constraint keys are documented in the `Network.constraint_info` when [`build_bfm!`](@ref) is
called. 


```julia
using BranchFlowModel
using CommonOPF
using JuMP
using ECOS


net = Network_Papavasiliou_2018()

m = JuMP.Model(ECOS.Optimizer) 

# modify the power balance constraints at bus 11, adding real and reactive power injection variables
b = "11"
@variable(m, 0.4 >= pgen11 >= 0)
@variable(m, 0.4 >= qgen11 >= 0)

JuMP.delete.(m, m[:power_balance_constraints][b][:real])

m[:power_balance_constraints][b][:real] = @constraint(m, 
    sum( m[:pij][(i,b)][1] for i in i_to_j(b, net) )
    - sum( m[:lij][(i,b)][1] * rij(i,b,net) for i in i_to_j(b, net) ) 
    + pgen11 - net[b][:Load].kws1[1] == 0
)

JuMP.delete.(m, m[:power_balance_constraints][b][:reactive])
m[:power_balance_constraints][b][:reactive] = @constraint(m, 
    sum( m[:pij][(i,b)][1] for i in i_to_j(b, net) )
    - sum( m[:lij][(i,b)][1] * rij(i,b,net) for i in i_to_j(b, net) ) 
    + qgen11 - net[b][:Load].kvars1[1] == 0
)
```

See the [JuMP documentation](https://jump.dev/JuMP.jl/stable/manual/constraints/#Delete-a-constraint) for more on deleting constraints.
