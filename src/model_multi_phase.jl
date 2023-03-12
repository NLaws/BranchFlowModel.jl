"""
    build_model!(m::JuMP.AbstractModel, p::Inputs{MultiPhase})

Add variables and constraints to `m` using the values in `p`. Calls the following functions:
```julia
add_variables(m, p)
constrain_power_balance(m, p)
constrain_substation_voltage(m, p)
constrain_KVL(m, p)
constrain_loads(m, p)
```
"""
function build_model!(m::JuMP.AbstractModel, p::Inputs{MultiPhase})
    add_variables(m, p)  # PSD done in add_variables
    constrain_power_balance(m, p)
    constrain_KVL(m, p)
    constrain_loads(m, p)
end

# hack for zero problem in MutableArithmetics
import Base: zero
zero(::Type{Any}) = 0.0
zero(::Type{Union{Float64, GenericAffExpr}}) = 0.0


"""
    add_variables(m, p::Inputs{MultiPhase})

Create complex variables:
- `m[:w]` are 3x3 Hermitian matrices of voltage squared (V*V^T)
- `m[:l]` are 3x3 Hermitian matrices of current squared (I*I^T)
- `m[:Sj]` are 3x1 matrices of net power injections (at bus j)
- `m[:Sij]` are 3x3 Complex matrices of line flow powers (from i to j)

The positive semi-definite constraints are also defined and stored as
- `m[:H][t][j]` where `t` is time step and `j` is the bus name.

All of the variable containers have typeof `Dict{Int, Dict{String, AbstractVecOrMat}}``.
- The first index is time step (integer)
- The second index is bus or line (string)
- and finally a matrix of complex variables


Some examples of using variables:
```julia

value.(m[:Sj][1]["671"])

value(variable_by_name(m, "real(Sj_1_671_1)"))

value.(m[:w][1]["671"])

value.(m[:Sj][1][p.substation_bus]) 

for b in keys(p.Pload)
    println(b, "  ", value.(m[:Sj][1][b]))
end

fix(variable_by_name(m, "real(Sj_1_645_3)"), 0.0, force=true)
```
"""
function add_variables(m, p::Inputs{MultiPhase})
    T = 1:p.Ntimesteps

    # type for inner dicts of variable containers,
    # which are dicts with time and bus keys
    S = Dict{String, AbstractVecOrMat}
    # voltage squared is Hermitian
    w = Dict{Int64, S}()
    # current squared is Hermitian
    l = Dict{Int64, S}()
    # complex line powers (at the sending end)
    Sij = Dict{Int64, S}()
    # complex powers injections 
    Sj = Dict{Int64, S}()

    m[:H] = Dict{Int64, S}()  # to store PSD matrices

    for t in T
        # fix head voltage at p.v0
        v0 = [p.v0 + 0im; -0.5*p.v0 + im*sqrt(3)/2*p.v0; -0.5*p.v0 - im*sqrt(3)/2*p.v0]
        w[t] = Dict(p.substation_bus => v0*cj(v0))
        # empty dicts for line values in each time step to fill
        l[t] = Dict()
        Sij[t] = Dict()
        # slack bus power injection
        Sj[t] = Dict(p.substation_bus =>  @variable(m, [1:3] in ComplexPlane(), 
            base_name="Sj_" * string(t) *"_"* p.substation_bus,
            lower_bound = p.P_lo_bound + p.Q_lo_bound*im, 
            upper_bound = p.P_up_bound + p.Q_up_bound*im
        ))
        m[:H][t] = Dict()
        
        for j in keys(p.phases_into_bus)
            if j == p.substation_bus
                # already filled Sj and w and there are no lines into substation_bus
                continue
            end

            i = i_to_j(j, p)[1]
            i_j = string(i * "-" * j)  # for radial network there is only one i in i_to_j

            # initialize line flows and injections as zeros that will remain for undefined phases

            # Sij is 3x3
            Sij[t][i_j] = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
            for phs1 in p.phases_into_bus[j], phs2 in p.phases_into_bus[j]
                Sij[t][i_j][phs1, phs2] = @variable(m, 
                    set = ComplexPlane(), base_name="Sij_" * string(t) *"_"* i_j *"_"*  string(phs1) * string(phs2), 
                    lower_bound = p.P_lo_bound + p.Q_lo_bound*im,
                    upper_bound = p.P_up_bound + p.Q_up_bound*im,
                )
            end

            # Sj is 3x1
            Sj[t][j] = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, reshape([0.0im; 0im; 0im], 3, 1))
            for phs in p.phases_into_bus[j]  # fill in Sj vector with variables for phases going into bus j
                real_load = 0.0
                if j in keys(p.Pload)
                    if phs in keys(p.Pload[j])
                        real_load = p.Pload[j][phs][t]
                    end
                end
        
                reactive_load = 0.0
                if j in keys(p.Qload)
                    if phs in keys(p.Qload[j])
                        reactive_load = p.Qload[j][phs][t]
                    end
                end
                Sj[t][j][phs] = @variable(m, 
                    set = ComplexPlane(), base_name="Sj_" * string(t) *"_"* j *"_"* string(phs), 
                    lower_bound = p.P_lo_bound + p.Q_lo_bound*im,
                    upper_bound = p.P_up_bound + p.Q_up_bound*im,
                    start = -real_load - reactive_load*im
                )
            end

            # Hermitian variables
            l[t][i_j] = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
            w[t][j]   = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])

            # fill in variables for all combinations of phases
            # TODO better bounds
            for phs1 in p.phases_into_bus[j], phs2 in p.phases_into_bus[j]

                l_up_bound = p.Isqaured_up_bounds[get_ijlinecode(i,j,p)]

                if phs1 <= phs2  # upper triangle of Hermitian matrices

                    if phs1 == phs2 # diagonal terms are real
                        # TODO THE BOUNDS LEAD TO INFEASIBLE PROBLEMS
                        l[t][i_j][phs1, phs2] = @variable(m, 
                            base_name="l_" * string(t) *"_"* i_j *"_"*  string(phs1) * string(phs2), 
                            # lower_bound = 0.0,  # diagonal value are real, positive
                            # upper_bound = l_up_bound
                        )

                        w[t][j][phs1, phs2] = @variable(m, 
                            base_name="w_" * string(t) *"_"* j *"_"* string(phs1) * string(phs2), 
                            lower_bound = p.v_lolim^2,
                            upper_bound = p.v_uplim^2,
                            start = 1.0
                        )

                    else
                        l[t][i_j][phs1, phs2] = @variable(m, 
                            set = ComplexPlane(), base_name="l_" * string(t) *"_"* i_j *"_"*  string(phs1) * string(phs2), 
                            # lower_bound = 0.0 + 0.0*im,  # must have negative imaginary parts in Hermitian matrix
                            # upper_bound = l_up_bound + l_up_bound*im
                        )

                        w[t][j][phs1, phs2] = @variable(m, 
                            set = ComplexPlane(), base_name="w_" * string(t) *"_"* j *"_"* string(phs1) * string(phs2), 
                            # lower_bound = p.v_lolim^2 + p.v_lolim^2*im,  # must have negative imaginary parts in Hermitian matrix
                            # upper_bound = p.v_uplim^2 + p.v_uplim^2*im,
                            start = 1.0 + 0.0im
                        )
                    end
                end
            end
            # TODO line limit inputs
            
            l[t][i_j] = Hermitian(l[t][i_j])
            w[t][j]   = Hermitian(w[t][j])

            M = Hermitian([
                w[t][j]          Sij[t][i_j];
                cj(Sij[t][i_j])  l[t][i_j]
            ])
            m[:H][t][j] = M

            @constraint(m, M in HermitianPSDCone())
        end
    end
    m[:w] = w
    m[:l] = l
    m[:Sij] = Sij
    m[:Sj] = Sj
    
    nothing
end


"""
    function constrain_power_balance(m, p::Inputs{MultiPhase})

Sij in - losses == sum of line flows out + net injection
NOTE: using sum over Pij for future expansion to mesh grids
i -> j -> k

All of the power balance constraints are stored in `m[:loadbalcons]` with the bus name (string)
as the first index. For example `m[:loadbalcons]["busname"]` will give the constrain container
from JuMP for all time steps.
"""
function constrain_power_balance(m, p::Inputs{MultiPhase})
    Sⱼ = m[:Sj]
    Sᵢⱼ = m[:Sij]
    lᵢⱼ = m[:l]
    m[:loadbalcons] = Dict()
    # TODO change Pⱼ and Qⱼ to expressions, make P₀ and Q₀ dv's, which will reduce # of variables
    # by (Nnodes - 1)*8760 and number of constraints by 6*(Nnodes - 1)*8760
    for j in p.busses
        if isempty(i_to_j(j, p)) && !isempty(j_to_k(j, p)) # source nodes, injection = flows out
            con = @constraint(m,  [t in 1:p.Ntimesteps],
                Sⱼ[t][j] .- diag( sum( Sᵢⱼ[t][string(j*"-"*k)] for k in j_to_k(j, p) ) ) .== 0
            )
        elseif isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # unconnected nodes
            @warn "Bus $j has no edges, setting Sⱼ and Qⱼ to zero."
            con = @constraint(m, [t in 1:p.Ntimesteps],
                diag(Sⱼ[t][j]) .== 0
            )
        elseif !isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # leaf nodes / sinks, flows in = draw out
            con = @constraint(m, [t in 1:p.Ntimesteps],
                diag(
                    sum( 
                        Sᵢⱼ[t][string(i*"-"*j)] .- zij(i,j,p) * lᵢⱼ[t][string(i*"-"*j)]
                    for i in i_to_j(j, p))
                )
                + Sⱼ[t][j] .== 0
            )
        else
            con =  @constraint(m, [t in 1:p.Ntimesteps],
                diag(
                    sum( 
                        Sᵢⱼ[t][string(i*"-"*j)] .- zij(i,j,p) * lᵢⱼ[t][string(i*"-"*j)]
                    for i in i_to_j(j, p))
                ) .+
                Sⱼ[t][j] .- 
                diag( sum( Sᵢⱼ[t][string(j*"-"*k)] for k in j_to_k(j, p) ) ) .== 0
            )
        end
        m[:loadbalcons][j] = con
    end
    # TODO Farivar and Low have b*v and q*v in these equations for shunts? neglected for now

    nothing
end


"""
    matrix_phases_to_vec(M::AbstractMatrix{T}, phases::AbstractVector{Int}) where T

Used in defining the KVL constraints, this method returns `M` for only
the indices in `phases`.
"""
function matrix_phases_to_vec(M::AbstractMatrix{T}, phases::AbstractVector{Int}) where T
    v = T[]
    for i in phases, j in phases 
        push!(v, M[i,j])
    end
    return v
end


"""
    constrain_KVL(m, p::Inputs{MultiPhase})

Add the voltage drop definintions between busses.
"""
function constrain_KVL(m, p::Inputs{MultiPhase})
    w = m[:w]
    Sᵢⱼ = m[:Sij]
    lᵢⱼ = m[:l]

    T = matrix_phases_to_vec  # "T" for transform
    m[:kvl] = Dict{String, AbstractArray}()
    for j in p.busses  # substation_bus in here but has empty i_to_j(j, p)
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            phs = p.phases_into_bus[j]
            i_j = string(i*"-"*j)
            z = zij(i,j,p)
            # need to slice w[t][i] by phases in i_j
            con = @constraint(m, [t in 1:p.Ntimesteps],
                T(w[t][j], phs) .== T(w[t][i], phs)
                    - T(Sᵢⱼ[t][i_j] * cj(z) + z * cj(Sᵢⱼ[t][i_j]), phs)
                    + T(z * lᵢⱼ[t][i_j] * cj(z), phs)
            );
            m[:kvl][j] = con
        end
    end
    nothing
end


"""
    constrain_loads(m, p::Inputs{MultiPhase})

- set loads to negative of Inputs.Pload and Inputs.Qload, 
    which are normalized by Sbase when creating Inputs.
- keys of Pload and Qload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus

Each of the power injection constraints are stored in the model under `m[:injectioncons]`.
To acces the constraints you can:
```julia
m[:injectioncons]["busname"][t,phs]
```
where `t` is the integer time step and `phs` is the integer phase.
"""
function constrain_loads(m, p::Inputs{BranchFlowModel.MultiPhase})
    Sⱼ = m[:Sj]
    m[:injectioncons] = Dict()
    for j in setdiff(p.busses, [p.substation_bus])

        # default to zero injection
        real_load = Dict(phs => zeros(p.Ntimesteps) for phs in [1,2,3])

        if j in keys(p.Pload)
            for phs in keys(p.Pload[j])
                real_load[phs] = p.Pload[j][phs]
            end
        end

        # default to zero injection
        reactive_load = Dict(phs => zeros(p.Ntimesteps) for phs in [1,2,3])

        if j in keys(p.Qload)
            for phs in keys(p.Qload[j])
                reactive_load[phs] = p.Qload[j][phs]
            end
        end

        con = @constraint(m, [t in 1:p.Ntimesteps, phs in p.phases_into_bus[j]],
            Sⱼ[t][j][phs] == (-real_load[phs][t] - reactive_load[phs][t]im ) / p.Sbase
        )
        m[:injectioncons][j] = con
    end
    nothing
end
