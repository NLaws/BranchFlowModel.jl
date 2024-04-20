"""
    build_model!(m::JuMP.AbstractModel, net::Network{MultiPhase})

Add variables and constraints to `m` using the values in `p`. Calls the following functions:
```julia
add_variables(m, p)
constrain_power_balance(m, p)
constrain_substation_voltage(m, p)
constrain_KVL(m, p)
constrain_loads(m, p)
```
"""
function build_model!(m::JuMP.AbstractModel, net::Network{MultiPhase})
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
    add_variables(m, net::Network{MultiPhase})

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
- and finally a matrix or vector of complex variables


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
function add_variables(m, net::Network{MultiPhase})
    # type for inner dicts of variable containers, which are dicts with time and bus keys
    S = Dict{String, AbstractVecOrMat}
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

    # fix head voltage at net.v0; if real values provided we assume 120deg phase shift
    if typeof(net.v0) <: Real
        v0 = [net.v0 + 0im; -0.5*net.v0 + im*sqrt(3)/2*net.v0; -0.5*net.v0 - im*sqrt(3)/2*net.v0]
        v0 =  v0*cj(v0)
    elseif typeof(net.v0) <: AbstractVector{<:Real}
        v0 = [
            net.v0[1] + 0im; 
            -0.5*net.v0[2] + im*sqrt(3)/2*net.v0[2]; 
            -0.5*net.v0[3] - im*sqrt(3)/2*net.v0[3]
        ]
        v0 = v0*cj(v0)
    elseif typeof(net.v0) <: AbstractVector{<:Complex}
        v0 = net.v0*cj(net.v0)
    else  # matrix provided
        v0 = net.v0
    end

    for t in 1:p.Ntimesteps
        m[:w][t] = Dict(p.substation_bus => v0)
        # empty dicts for line values in each time step to fill
        m[:l][t] = Dict()
        m[:Sij][t] = Dict()
        # slack bus power injection
        m[:Sj][t] = Dict(p.substation_bus =>  @variable(m, [1:3] in ComplexPlane(), 
            base_name="Sj_" * string(t) *"_"* p.substation_bus,
            # lower_bound = p.P_lo_bound + p.Q_lo_bound*im, 
            # upper_bound = p.P_up_bound + p.Q_up_bound*im
        ))
        m[:H][t] = Dict()

        # inner method to loop over
        function define_vars_downstream(i::String, t::Int, m::JuMP.AbstractModel, net::Network)
            for j in j_to_k(i, net)  # i -> j -> k
                i_j = string(i * "-" * j)  # for radial network there is only one i in i_to_j

                # initialize line flows and injections as zeros that will remain for undefined phases
                # Sij is 3x3
                m[:Sij][t][i_j] = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
                for phs1 in p.phases_into_bus[j], phs2 in p.phases_into_bus[j]
                    m[:Sij][t][i_j][phs1, phs2] = @variable(m, 
                        set = ComplexPlane(), base_name="Sij_" * string(t) *"_"* i_j *"_"*  string(phs1) * string(phs2), 
                        # lower_bound = p.P_lo_bound + p.Q_lo_bound*im,
                        # upper_bound = p.P_up_bound + p.Q_up_bound*im,
                    )
                end

                # Sj is 3x1
                m[:Sj][t][j] = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, reshape([0.0im; 0im; 0im], 3, 1))
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
                    m[:Sj][t][j][phs] = @variable(m, 
                        set = ComplexPlane(), base_name="Sj_" * string(t) *"_"* j *"_"* string(phs), 
                        lower_bound = p.P_lo_bound + p.Q_lo_bound*im,
                        upper_bound = p.P_up_bound + p.Q_up_bound*im,
                        start = -real_load - reactive_load*im
                    )
                end

                # Hermitian variables
                m[:l][t][i_j] = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
                m[:w][t][j]   = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])

                # fill in variables for all combinations of phases
                # TODO better bounds
                for phs1 in p.phases_into_bus[j], phs2 in p.phases_into_bus[j]

                    l_up_bound = p.Isquared_up_bounds[get_ijlinecode(i,j,p)]

                    if phs1 <= phs2  # upper triangle of Hermitian matrices

                        if phs1 == phs2 # diagonal terms are real
                            # TODO THE BOUNDS LEAD TO INFEASIBLE PROBLEMS
                            m[:l][t][i_j][phs1, phs2] = @variable(m, 
                                base_name="l_" * string(t) *"_"* i_j *"_"*  string(phs1) * string(phs2), 
                                # lower_bound = 0.0,  # diagonal value are real, positive
                                # upper_bound = l_up_bound
                            )

                            m[:w][t][j][phs1, phs2] = @variable(m, 
                                base_name="w_" * string(t) *"_"* j *"_"* string(phs1) * string(phs2), 
                                lower_bound = p.v_lolim^2,
                                upper_bound = p.v_uplim^2,
                                start = 1.0
                            )

                        else
                            m[:l][t][i_j][phs1, phs2] = @variable(m, 
                                set = ComplexPlane(), base_name="l_" * string(t) *"_"* i_j *"_"*  string(phs1) * string(phs2), 
                                # lower_bound = 0.0 + 0.0*im,  # must have negative imaginary parts in Hermitian matrix
                                # upper_bound = l_up_bound + l_up_bound*im
                            )

                            m[:w][t][j][phs1, phs2] = @variable(m, 
                                set = ComplexPlane(), base_name="w_" * string(t) *"_"* j *"_"* string(phs1) * string(phs2), 
                                # lower_bound = p.v_lolim^2 + p.v_lolim^2*im,  # must have negative imaginary parts in Hermitian matrix
                                # upper_bound = p.v_uplim^2 + p.v_uplim^2*im,
                                start = 1.0 + 0.0im
                            )
                        end
                    end
                end
                # TODO line limit inputs
                
                m[:l][t][i_j] = Hermitian(m[:l][t][i_j])
                m[:w][t][j]   = Hermitian(m[:w][t][j])

                w_i =  phi_ij(j, p, m[:w][t][i])

                M = Hermitian([
                    w_i                  m[:Sij][t][i_j];
                    cj(m[:Sij][t][i_j])  m[:l][t][i_j]
                ])
                m[:H][t][j] = M

                @constraint(m, M in HermitianPSDCone())
            end
            nothing
        end
        
        # have to traverse down the tree in order b/c the semi-definite constraints have i and i_j variables
        i = p.substation_bus
        define_vars_downstream(i, t, m, p)

        function recursive_variables(i::String, t::Int, m::JuMP.AbstractModel, net::Network)
            for j in j_to_k(i, net)
                define_vars_downstream(j, t, m, p)
                recursive_variables(j, t, m, p)
            end
        end
        recursive_variables(i, t, m, p)
    end
    
    nothing
end


"""
    function constrain_power_balance(m, net::Network{MultiPhase})

Sij in - losses == sum of line flows out + net injection
NOTE: using sum over Pij for future expansion to mesh grids
i -> j -> k

All of the power balance constraints are stored in `m[:loadbalcons]` with the bus name (string)
as the first index. For example `m[:loadbalcons]["busname"]` will give the constrain container
from JuMP for all time steps.
"""
function constrain_power_balance(m, net::Network{MultiPhase})
    Sj = m[:Sj]
    Sij = m[:Sij]
    lij = m[:l]
    m[:loadbalcons] = Dict()
    # TODO change Pj and Qj to expressions, make P₀ and Q₀ dv's, which will reduce # of variables
    # by (Nnodes - 1)*8760 and number of constraints by 6*(Nnodes - 1)*8760
    for j in p.busses
        if isempty(i_to_j(j, p)) && !isempty(j_to_k(j, p)) # source nodes, injection = flows out
            con = @constraint(m,  [t in 1:p.Ntimesteps],
                Sj[t][j] - sum( diag( Sij[t][string(j*"-"*k)] ) for k in j_to_k(j, p) ) .== 0
            )
        elseif isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # unconnected nodes
            @warn "Bus $j has no edges, setting Sj and Qj to zero."
            con = @constraint(m, [t in 1:p.Ntimesteps],
                diag(Sj[t][j]) .== 0
            )
        elseif !isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # leaf nodes / sinks, flows in = draw out
            con = @constraint(m, [t in 1:p.Ntimesteps],
                sum( diag( 
                    Sij[t][string(i*"-"*j)] - zij(i,j,p) * lij[t][string(i*"-"*j)]
                ) for i in i_to_j(j, p) )
                + Sj[t][j] .== 0
            )
        else  # node with lines in and out
            con =  @constraint(m, [t in 1:p.Ntimesteps],
                sum( diag( 
                    Sij[t][string(i*"-"*j)] - zij(i,j,p) * lij[t][string(i*"-"*j)]
                ) for i in i_to_j(j, p) )
                + Sj[t][j]
                - sum( diag( Sij[t][string(j*"-"*k)] ) for k in j_to_k(j, p) ) .== 0
            )
        end
        m[:loadbalcons][j] = con
    end
    # TODO add shunts, diag( openDSS-cmatrix * m[:w][t][j] ) ?

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
    constrain_KVL(m, net::Network{MultiPhase})

Add the voltage drop definintions between busses.
"""
function constrain_KVL(m, net::Network{MultiPhase})
    w = m[:w]
    Sij = m[:Sij]
    lij = m[:l]

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
                    - T(Sij[t][i_j] * cj(z) + z * cj(Sij[t][i_j]), phs)
                    + T(z * lij[t][i_j] * cj(z), phs)
            );
            m[:kvl][j] = con
        end
    end
    nothing
end


"""
    constrain_loads(m, net::Network{MultiPhase})

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
function constrain_loads(m, net::Network{BranchFlowModel.MultiPhase})
    Sj = m[:Sj]
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
            Sj[t][j][phs] == (-real_load[phs][t] - reactive_load[phs][t]im ) / p.Sbase
        )
        m[:injectioncons][j] = con
    end
    nothing
end
