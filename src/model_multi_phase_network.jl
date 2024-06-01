"""
    build_model!(m::JuMP.AbstractModel, net::Network{MultiPhase})

Add variables and constraints to `m` using the values in `net`. Calls the following functions:
```julia
add_variables(m, net)
constrain_power_balance(m, net)
constrain_substation_voltage(m, net)
constrain_KVL(m, net)
```
"""
function build_model!(m::JuMP.AbstractModel, net::Network{MultiPhase})
    add_variables(m, net)  # PSD done in add_variables
    constrain_power_balance(m, net)
    constrain_KVL(m, net)
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

value.(m[:Sj][1][net.substation_bus]) 

for b in real_load_busses(net)
    println(b, "  ", value.(m[:Sj][1][b]))
end

fix(variable_by_name(m, "real(Sj_1_645_3)"), 0.0, force=true)
```
"""
function add_variables(m, net::Network{MultiPhase})
    # TODO complex model containers in CommonOPF
    # type for inner dicts of variable containers, which are dicts with time and bus keys
    S_bus = Dict{String, AbstractVecOrMat}
    S_edge = Dict{Tuple{String, String}, AbstractVecOrMat}
    # voltage squared is Hermitian
    m[:w] = Dict{Int64, S_bus}()
    # current squared is Hermitian
    m[:l] = Dict{Int64, S_edge}()
    # complex line powers (at the sending end)
    m[:Sij] = Dict{Int64, S_edge}()
    # complex net powers injections 
    m[:Sj] = Dict{Int64, S_bus}()
    # Hermitian PSD matrices
    m[:H] = Dict{Int64, S_bus}()

    # fix head voltage at net.v0; if real values provided we assume 120deg phase shift
    if typeof(net.v0) <: Real
        v0 = [net.v0 + 0im; -0.5*net.v0 + im*sqrt(3)/2 * net.v0; -0.5*net.v0 - im*sqrt(3)/2 * net.v0]
        v0 =  v0*cj(v0)
    elseif typeof(net.v0) <: AbstractVector{<:Real}
        v0 = [
            net.v0[1] + 0im; 
            -0.5 * net.v0[2] + im*sqrt(3)/2 * net.v0[2]; 
            -0.5 * net.v0[3] - im*sqrt(3)/2 * net.v0[3]
        ]
        v0 = v0*cj(v0)
    elseif typeof(net.v0) <: AbstractVector{<:Complex}
        v0 = net.v0 * cj(net.v0)
    else  # matrix provided
        v0 = net.v0
    end

    for t in 1:net.Ntimesteps
        m[:w][t] = Dict(net.substation_bus => v0 * cj(v0))
        # empty dicts for line values in each time step to fill
        m[:l][t] = Dict()
        m[:Sij][t] = Dict()
        # slack bus power injection
        m[:Sj][t] = Dict(net.substation_bus =>  @variable(m, [1:3] in ComplexPlane(), 
            base_name="Sj_" * string(t) *"_"* net.substation_bus,
            # lower_bound = p.P_lo_bound + p.Q_lo_bound*im, 
            # upper_bound = p.P_up_bound + p.Q_up_bound*im
        ))
        m[:H][t] = Dict()

        # inner method to loop over
        function define_vars_downstream(i::String, t::Int, m::JuMP.AbstractModel, net::Network)
            for j in j_to_k(i, net)  # i -> j -> k
                i_j = (i, j)  # for radial network there is only one i in i_to_j

                # initialize line flows and injections as zeros that will remain for undefined phases
                # Sij is 3x3
                m[:Sij][t][i_j] = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
                for phs1 in phases_into_bus(net, j), phs2 in phases_into_bus(net, j)
                    m[:Sij][t][i_j][phs1, phs2] = @variable(m, 
                        set = ComplexPlane(), base_name="Sij_" * string(t) *"_"* string(i) *"_"* string(j) *"_"*  string(phs1) * string(phs2), 
                        # lower_bound = p.P_lo_bound + p.Q_lo_bound*im,
                        # upper_bound = p.P_up_bound + p.Q_up_bound*im,
                    )
                end

                # Hermitian variables
                m[:l][t][i_j] = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
                m[:w][t][j]   = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])

                # fill in variables for all combinations of phases
                # TODO better bounds
                for phs1 in phases_into_bus(net, j), phs2 in phases_into_bus(net, j)

                    # l_up_bound = p.Isquared_up_bounds[get_ijlinecode(i,j,p)]

                    if phs1 <= phs2  # upper triangle of Hermitian matrices

                        if phs1 == phs2 # diagonal terms are real
                            # TODO THE BOUNDS LEAD TO INFEASIBLE PROBLEMS
                            m[:l][t][i_j][phs1, phs2] = @variable(m, 
                                base_name="l_" * string(t) *"_"* string(i) *"_"* string(j) *"_"*  string(phs1) * string(phs2), 
                                # lower_bound = 0.0,  # diagonal value are real, positive
                                # upper_bound = l_up_bound
                            )

                            m[:w][t][j][phs1, phs2] = @variable(m, 
                                base_name="w_" * string(t) *"_"* j *"_"* string(phs1) * string(phs2), 
                                lower_bound = net.v_lolim^2,
                                upper_bound = net.v_uplim^2,
                                start = 1.0
                            )

                        else
                            m[:l][t][i_j][phs1, phs2] = @variable(m, 
                                set = ComplexPlane(), base_name="l_" * string(t) *"_"* string(i) *"_"* string(j) *"_"*  string(phs1) * string(phs2), 
                                # lower_bound = 0.0 + 0.0*im,  # must have negative imaginary parts in Hermitian matrix
                                # upper_bound = l_up_bound + l_up_bound*im
                            )

                            m[:w][t][j][phs1, phs2] = @variable(m, 
                                set = ComplexPlane(), base_name="w_" * string(t) *"_"* j *"_"* string(phs1) * string(phs2), 
                                # lower_bound = net.v_lolim^2 + net.v_lolim^2*im,  # must have negative imaginary parts in Hermitian matrix
                                # upper_bound = net.v_uplim^2 + net.v_uplim^2*im,
                                start = 1.0 + 0.0im
                            )
                        end
                    end
                end
                # TODO line limit inputs
                
                m[:l][t][i_j] = Hermitian(m[:l][t][i_j])
                m[:w][t][j]   = Hermitian(m[:w][t][j])

                w_i =  phi_ij(j, net, m[:w][t][i])

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
        i = net.substation_bus
        define_vars_downstream(i, t, m, net)

        function recursive_variables(j::String, t::Int, m::JuMP.AbstractModel, net::Network)
            for k in j_to_k(j, net)
                define_vars_downstream(k, t, m, net)
                recursive_variables(k, t, m, net)
            end
        end
        recursive_variables(i, t, m, net)
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
    Sij = m[:Sij]
    lij = m[:l]
    m[:loadbalcons] = Dict()
    
    for j in busses(net)

        # check for loads
        Pj = [zeros(net.Ntimesteps) for _ in 1:3] # first dim is phase, like Pj[phs][t]
        Qj = [zeros(net.Ntimesteps) for _ in 1:3]
        if j in real_load_busses(net)
            for phs in 1:3  # put an Sj method in CommonOPF? to make complex vector of phases
                Pj[phs] = -net[j, :kws, phs] * 1e3 / net.Sbase
            end
        end
        if j in reactive_load_busses(net)
            for phs in 1:3 
                Qj[phs] = -net[j, :kvars, phs] * 1e3 / net.Sbase
            end
        end
        Sj = Pj + im * Qj
        Sj = hcat(Sj...)  # time X phase

        # source nodes, injection = flows out
        if isempty(i_to_j(j, net)) && !isempty(j_to_k(j, net))

            if j == net.substation_bus   # include the slack power variables
                m[:loadbalcons][j] = @constraint(m,  [t in 1:net.Ntimesteps],
                    m[:Sj][t][j] + Sj[t, :] - sum( diag( Sij[t][(j,k)] ) for k in j_to_k(j, net) ) .== 0
                )
            else  # a source node with known injection
                m[:loadbalcons][j] = @constraint(m,  [t in 1:net.Ntimesteps],
                    Sj[t, :] - sum( diag( Sij[t][(j,k)] ) for k in j_to_k(j, net) ) .== 0
                )

            end
        
        # unconnected nodes
        elseif isempty(i_to_j(j, net)) && isempty(j_to_k(j, net))
            @warn "Bus $j has no edges. Setting loadbalcons to zeros"
            m[:loadbalcons][j] = zeros(net.Ntimesteps)

        # leaf nodes / sinks, flows in = draw out
        elseif !isempty(i_to_j(j, net)) && isempty(j_to_k(j, net))
            m[:loadbalcons][j] = @constraint(m, [t in 1:net.Ntimesteps],
                sum( diag( 
                    Sij[t][(i,j)] - zij_per_unit(i,j,net) * lij[t][(i,j)]
                ) for i in i_to_j(j, net) )
                + Sj[t, :] .== 0
            )

        # node with lines in and out
        else
            m[:loadbalcons][j] = @constraint(m, [t in 1:net.Ntimesteps],
                sum( diag( 
                    Sij[t][(i,j)] - zij_per_unit(i,j,net) * lij[t][(i,j)]
                ) for i in i_to_j(j, net) )
                + Sj[t, :]
                - sum( diag( Sij[t][(j,k)] ) for k in j_to_k(j, net) ) .== 0
            )
        end
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
    for j in busses(net)  # substation_bus in here but has empty i_to_j(j, net)
        for i in i_to_j(j, net)  # for radial network there is only one i in i_to_j
            phases = phases_into_bus(net, j)
            z = zij_per_unit(i,j,net)
            # need to slice w[t][i] by phases in edge (i, j)
            con = @constraint(m, [t in 1:net.Ntimesteps],
                T( w[t][j], phases ) .== T( w[t][i], phases )
                    - T(
                        Sij[t][(i, j)] * cj(z) + z * cj(Sij[t][(i, j)]),
                    phases)
                    + T(
                        z * lij[t][(i, j)] * cj(z), 
                    phases)
            );
            m[:kvl][j] = con
        end
    end
    nothing
end
