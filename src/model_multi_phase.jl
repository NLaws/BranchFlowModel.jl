"""
    build_bfm!(m::JuMP.AbstractModel, net::Network{MultiPhase}, mtype::ModelType=Unrelaxed)

Top-level multiphase builder that dispatches the ModelType enum
"""
function build_bfm!(m::JuMP.AbstractModel, net::Network{MultiPhase}, mtype::ModelType=Unrelaxed)
    build_bfm!(m::JuMP.AbstractModel, net::Network{MultiPhase}, Val(mtype))
end


"""
    build_bfm!(m::JuMP.AbstractModel, net::Network{MultiPhase}, ::Val{Unrelaxed})

Add variables and constraints to `m` using the values in `net` to make an unrelaxed branch flow
model. Calls the following functions:
```julia
add_bfm_variables(m, net)
constrain_bfm_nlp(m, net)
```
"""
function build_bfm!(m::JuMP.AbstractModel, net::Network{MultiPhase}, ::Val{Unrelaxed})
    add_bfm_variables(m, net)
    constrain_bfm_nlp(m, net)
end


"""
    build_bfm!(m::JuMP.AbstractModel, net::Network{MultiPhase}, ::Val{Semidefinite})

Add variables and constraints to `m` using the values in `net`. Calls the following functions:
```julia
add_sdp_variables(m, net)
constrain_power_balance(m, net)
constrain_KVL(m, net)
```
"""
function build_bfm!(m::JuMP.AbstractModel, net::Network{MultiPhase}, ::Val{Semidefinite})
    add_sdp_variables(m, net)  # PSD constraints done in add_sdp_variables
    constrain_power_balance(m, net)
    constrain_KVL(m, net)
end


# hack for zero problem in MutableArithmetics
import Base: zero
zero(::Type{Any}) = 0.0
zero(::Type{Union{Float64, GenericAffExpr}}) = 0.0


function substation_voltage(net::Network{MultiPhase})::Vector{ComplexF64}

    if typeof(net.v0) <: Real
        return [
            net.v0 + 0im; 
            -0.5*net.v0 - im*sqrt(3)/2 * net.v0; 
            -0.5*net.v0 + im*sqrt(3)/2 * net.v0
        ]

    elseif typeof(net.v0) <: AbstractVector{<:Real}
        return [
            net.v0[1] + 0im; 
            -0.5 * net.v0[2] - im*sqrt(3)/2 * net.v0[2]; 
            -0.5 * net.v0[3] + im*sqrt(3)/2 * net.v0[3]
        ]

    elseif typeof(net.v0) <: AbstractVector{<:Complex}
        return net.v0

    else  
        throw(@error "unsupported type for Network.v0 $(typeof(net.v0))")
    end

    return w0
end


"""
    function substation_voltage_squared(net::Network{MultiPhase})::Matrix{ComplexF64}

Consruct voltage matrix for substation from `net.v0`. If a real value or vector of real values is
provided then a 120 degree phase shift is assumed. A vector of complex values is also supported.

"""
function substation_voltage_squared(net::Network{MultiPhase})::Matrix{ComplexF64}

    if typeof(net.v0) <: Real
        v0 = [net.v0 + 0im; -0.5*net.v0 + im*sqrt(3)/2 * net.v0; -0.5*net.v0 - im*sqrt(3)/2 * net.v0]
        w0 =  v0 * cj(v0)

    elseif typeof(net.v0) <: AbstractVector{<:Real}
        v0 = [
            net.v0[1] + 0im; 
            -0.5 * net.v0[2] + im*sqrt(3)/2 * net.v0[2]; 
            -0.5 * net.v0[3] - im*sqrt(3)/2 * net.v0[3]
        ]
        w0 = v0 * cj(v0)

    elseif typeof(net.v0) <: AbstractVector{<:Complex}
        w0 = net.v0 * cj(net.v0)

    else  
        throw(@error "unsupported type for Network.v0 $(typeof(net.v0))")
    end

    return w0
end


"""
    add_sdp_variables(m, net::Network{MultiPhase})

Create complex variables:
- `m[:w]` are 3x3 Hermitian matrices of voltage squared (V*V^T)
- `m[:l]` are 3x3 Hermitian matrices of current squared (I*I^T)
- `m[:Sj]` are 3x1 matrices of net power injections (at bus j)
- `m[:Sij]` are 3x3 Complex matrices of line flow powers (from i to j)

If PSD is `true` then the positive semi-definite constraints are also defined and stored as
- `m[:H][j][t]` where `t` is time step and `j` is the bus name.

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
function add_sdp_variables(m, net::Network{MultiPhase})
    # type for inner dicts of variable containers, which are dicts with bus and time keys
    # voltage squared is Hermitian
    m[:w] = multiphase_bus_variable_container()
    # current squared is Hermitian
    m[:l] = multiphase_edge_variable_container()
    # complex line powers (at the sending end)
    m[:Sij] = multiphase_edge_variable_container()
    # complex net powers injections 
    m[:Sj] = multiphase_bus_variable_container()
    # Hermitian PSD matrices
    m[:H] = multiphase_bus_variable_container()

    for t in 1:net.Ntimesteps
        m[:w][net.substation_bus][t] = substation_voltage_squared(net)
        # slack bus power injection
        m[:Sj][net.substation_bus][t] = @variable(m, [1:3] in ComplexPlane(), 
            base_name="Sj_" * string(t) *"_"* net.substation_bus,
        )
        if !ismissing(net.bounds.s_lower_real)
            @constraint(m, net.bounds.s_lower_real .<= real( m[:Sj][net.substation_bus][t] ))
        end
        if !ismissing(net.bounds.s_upper_real)
            @constraint(m, real( m[:Sj][net.substation_bus][t] ) .<= net.bounds.s_upper_real)
        end
        if !ismissing(net.bounds.s_lower_imag)
            @constraint(m, net.bounds.s_lower_imag .<= imag( m[:Sj][net.substation_bus][t] ))
        end
        if !ismissing(net.bounds.s_upper_imag)
            @constraint(m, imag( m[:Sj][net.substation_bus][t] ) .<= net.bounds.s_upper_imag)
        end

        # inner method to loop over
        function define_vars_downstream(i::String, t::Int, m::JuMP.AbstractModel, net::Network)
            for j in j_to_k(i, net)  # i -> j -> k
                i_j = (i, j)  # for radial network there is only one i in i_to_j

                # initialize line flows and injections as zeros that will remain for undefined phases
                # Sij is 3x3
                m[:Sij][i_j][t] = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
                for phs1 in phases_into_bus(net, j), phs2 in phases_into_bus(net, j)
                    m[:Sij][i_j][t][phs1, phs2] = @variable(m, 
                        set = ComplexPlane(), base_name="Sij_" * string(t) *"_"* string(i) *"_"* string(j) *"_"*  string(phs1) * string(phs2), 
                        # upper_bound = net.bounds.s_upper + net.bounds.s_upper*im, 
                        # lower_bound = net.bounds.s_lower + net.bounds.s_lower*im
                    )
                end

                # Hermitian variables
                m[:l][i_j][t] = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
                m[:w][j][t]   = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])

                # fill in variables for all combinations of phases
                for phs1 in phases_into_bus(net, j), phs2 in phases_into_bus(net, j)

                    if phs1 <= phs2  # upper triangle of Hermitian matrices

                        # real diagonal terms
                        if phs1 == phs2
                            m[:l][i_j][t][phs1, phs2] = @variable(m, 
                                base_name="l_" * string(t) *"_"* string(i) *"_"* string(j) *"_"*  string(phs1) * string(phs2),
                                # upper_bound = net.bounds.i_upper_real^2,  # TODO bounds can be missing
                                lower_bound = 0.0,  # diagonal values are real, positive
                                start = 0.01
                            )

                            m[:w][j][t][phs1, phs2] = @variable(m, 
                                base_name="w_" * string(t) *"_"* j *"_"* string(phs1) * string(phs2),
                                # upper_bound = net.bounds.v_upper^2,
                                # lower_bound = net.bounds.v_lower^2,
                                start = 1.0
                            )
                        # complex off-diagonal terms
                        else
                            m[:l][i_j][t][phs1, phs2] = @variable(m, 
                                set = ComplexPlane(), base_name="l_" * string(t) *"_"* string(i) *"_"* string(j) *"_"*  string(phs1) * string(phs2),
                                # upper_bound =  net.bounds.i_upper^2 + im * net.bounds.i_upper^2,
                                # lower_bound =  net.bounds.i_lower^2 + im * net.bounds.i_lower^2,  # must have negative imaginary parts in Hermitian matrix
                            )

                            m[:w][j][t][phs1, phs2] = @variable(m, 
                                set = ComplexPlane(), base_name="w_" * string(t) *"_"* j *"_"* string(phs1) * string(phs2),
                                # upper_bound = net.bounds.v_upper^2 + net.bounds.v_upper^2*im,
                                # lower_bound = net.bounds.v_lower^2 + net.bounds.v_lower^2*im,  # must have negative imaginary parts in Hermitian matrix
                                start = 1.0 + 0.0im
                            )
                        end
                    end
                end
                # TODO line limit inputs
                
                m[:l][i_j][t] = Hermitian(m[:l][i_j][t])
                m[:w][j][t]   = Hermitian(m[:w][j][t])

                w_i =  phi_ij(j, net, m[:w][i][t])

                M = Hermitian([
                    w_i                  m[:Sij][i_j][t];
                    cj(m[:Sij][i_j][t])  m[:l][i_j][t]
                ])

                m[:H][j][t] = M

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
    add_bfm_variables(m, net::Network{MultiPhase})

Define complex variables for:
- `:v` bus voltages
- `:i` branch currents
- `:sij` branch power flows
- `:sj` bus net power injections

"""
function add_bfm_variables(m, net::Network{MultiPhase})
    # type for inner dicts of variable containers, which are dicts with time and bus keys
    # bus voltage vectors by time and phase
    m[:v] = multiphase_bus_variable_container()
    # sending current vectors by time and phase
    m[:i] = multiphase_edge_variable_container()
    # sending line power matrices = v_i i_ij^H
    m[:Sij] = multiphase_edge_variable_container()
    # bus net power injection vectors
    m[:Sj] = multiphase_bus_variable_container()

    for t in 1:net.Ntimesteps
        m[:v][net.substation_bus][t] = substation_voltage(net)
        # slack bus power injection
        m[:Sj][net.substation_bus][t] = @variable(m, [1:3] in ComplexPlane(), 
            base_name="Sj_" * string(t) *"_"* net.substation_bus,
            # upper_bound = net.bounds.s_upper + net.bounds.s_upper*im, 
            # lower_bound = net.bounds.s_lower + net.bounds.s_lower*im
        )

        # inner method to loop over
        function define_vars_downstream(i::String, t::Int, m::JuMP.AbstractModel, net::Network)
            for j in j_to_k(i, net)  # i -> j -> k
                i_j = (i, j)  # for radial network there is only one i in i_to_j

                # initialize line flows and injections as zeros that will remain for undefined phases
                # Sij is 3x3
                m[:Sij][i_j][t] = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])

                for phs1 in phases_into_bus(net, j), phs2 in phases_into_bus(net, j)
                    m[:Sij][i_j][t][phs1, phs2] = @variable(m, 
                        set = ComplexPlane(), 
                        # base_name is Sij_t_i_j_phs1phs2 like Sij_1_bus1_bus2_23
                        base_name="Sij_" * string(t) *"_"* string(i) *"_"* string(j) *"_"*  string(phs1) * string(phs2), 
                    )
                end

                # 3-element vectors for current and voltage
                add_complex_vector_of_phase_variable!(m, net, i_j, :i, t;
                    upper_bound_mag = net.bounds.i_upper_mag,
                    lower_bound_mag = net.bounds.i_lower_mag,
                )
                add_complex_vector_of_phase_variable!(m, net, j, :v, t;
                    upper_bound_mag = net.bounds.v_upper_mag,
                    lower_bound_mag = net.bounds.v_lower_mag,
                )
            end
            nothing
        end
        
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
    add_linear_variables(m, net::Network{MultiPhase})

Define real phase-vector variables for:
- `:vsqrd` bus voltage magnitude squared
- `:pij`, `:qij` branch power flows
- `:p0`, `:q0` the slack bus net power injection
"""
function add_linear_variables(m, net::Network{MultiPhase})
    es = edges(net)


    d = p.phases_into_bus
    d[p.substation_bus] = [1,2,3]
    T = 1:p.Ntimesteps

    # # bus injections
    # @variables m begin
    #     p.P_lo_bound <= Pj[b in p.busses, d[b], T] <= p.P_up_bound
    #     p.Q_lo_bound <= Qj[b in p.busses, d[b], T] <= p.Q_up_bound
    # end
    
    # voltage squared
    @variable(m, p.v_lolim^2 <= vsqrd[b in p.busses, d[b], T] <= p.v_uplim^2 ) 
    
    edge2phases = Dict(k=>v for (k,v) in zip(es, p.phases))  # this will fail for mesh network

    # line flows, power sent from i to j
    @variable(m, p.P_lo_bound <= Pij[e in es, edge2phases[e], T] <= p.P_up_bound )
    
    @variable(m, p.Q_lo_bound <= Qij[e in es, edge2phases[e], T] <= p.Q_up_bound )
    # TODO line flow limit inputs
    nothing
end


"""
    constrain_bfm_nlp(m, net::Network{MultiPhase})

Add the unrelaxed branch flow constraints.
"""
function constrain_bfm_nlp(m, net::Network{MultiPhase})
    v = m[:v]
    i_ij = m[:i]
    S_ij = m[:Sij]
    for t in 1:net.Ntimesteps

        for j in busses(net)
            for i in i_to_j(j, net)
            
                # sending end power
                @constraint(m, [t in 1:net.Ntimesteps],
                    S_ij[(i,j)][t] .== ( 
                        phi_ij(j, net, v[i][t]) * cj(i_ij[(i,j)][t]) 
                    )
                )

                # KVL
                @constraint(m, [t in 1:net.Ntimesteps],
                    phi_ij(j, net, v[i][t]) - v[j][t] .== zij_per_unit(i, j, net) * i_ij[(i,j)][t]
                )

            end

        end

    end
    constrain_power_balance(m, net)
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

    if :l in keys(m.obj_dict)
        Lij = m[:l]
    else
        Lij = multiphase_edge_variable_container()
        for t in 1:net.Ntimesteps
            for (i,j) in edges(net)
                Lij[(i,j)][t] = m[:i][(i,j)][t] * cj(m[:i][(i,j)][t])
            end
        end
    end

    if :w in keys(m.obj_dict)
        w = m[:w]
    else
        w = multiphase_bus_variable_container()
        for t in 1:net.Ntimesteps
            for j in busses(net)
                w[j][t] = m[:v][j][t] * cj(m[:v][j][t])
            end
        end
    end
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
                    m[:Sj][j][t] + Sj[t, :]
                    # - diag(w[j][t] * conj(yj(j, net))) * net.Zbase  # put yj in per-unit
                    - sum( diag( Sij[(j,k)][t] ) for k in j_to_k(j, net) )
                    .== 0
                )

            else  # a source node with known injection
                m[:loadbalcons][j] = @constraint(m,  [t in 1:net.Ntimesteps],
                    Sj[t, :] 
                    # - diag(w[j][t] * conj(yj(j, net))) * net.Zbase  # put yj in per-unit
                    - sum( diag( Sij[(j,k)][t] ) for k in j_to_k(j, net) ) 
                    .== 0
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
                    Sij[(i,j)][t] - zij_per_unit(i,j,net) * Lij[(i,j)][t]
                ) for i in i_to_j(j, net) )
                + Sj[t, :] 
                # - diag(w[j][t] * conj(yj(j, net))) * net.Zbase  # put yj in per-unit
                .== 0
            )

        # node with lines in and out
        else
            m[:loadbalcons][j] = @constraint(m, [t in 1:net.Ntimesteps],
                sum( diag( 
                    Sij[(i,j)][t] - zij_per_unit(i,j,net) * Lij[(i,j)][t]
                ) for i in i_to_j(j, net) )
                + Sj[t, :]
                # - diag(w[j][t] * conj(yj(j, net))) * net.Zbase  # put yj in per-unit
                - sum( diag( Sij[(j,k)][t] ) for k in j_to_k(j, net) ) 
                .== 0
            )
        end
    end

    nothing
end


"""
    matrix_phases_to_vec(M::AbstractMatrix{T}, phases::AbstractVector{Int}) where T

Used in defining the KVL constraints, this method returns the entries of `M` at the indices in 
`phases` in a vector.
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

Add the voltage drop definitions between busses.

``
    w_j = w_i - S_{ij} Z^{\\star} - Z S_{ij}^{\\star} + Z L_{ij} Z^{\\star}
``
"""
function constrain_KVL(m, net::Network{MultiPhase})
    w = m[:w]
    Sij = m[:Sij]
    Lij = m[:l]

    T = matrix_phases_to_vec  # "T" for transform
    m[:kvl] = Dict{String, AbstractArray}()
    for j in busses(net)  # substation_bus in here but has empty i_to_j(j, net)

        for i in i_to_j(j, net)  # for radial network there is only one i in i_to_j
            phases = phases_into_bus(net, j)

            if !( isa(net[(i,j)], CommonOPF.VoltageRegulator) )
                Z = zij_per_unit(i,j,net)
                # slice w[i][t] by phases in edge (i, j)
                m[:kvl][j] = @constraint(m, [t in 1:net.Ntimesteps],
                    T( w[j][t], phases ) .== T( w[i][t]
                        - (Sij[(i, j)][t] * cj(Z)
                           + Z * cj(Sij[(i, j)][t]))
                        + Z * Lij[(i, j)][t] * cj(Z), phases )
                );

            else
                # TODO account for phase angles/shifts
                reg = net[(i,j)]
                if !ismissing(reg.vreg_pu)
                    m[:kvl][j] = @constraint(m, [t in 1:net.Ntimesteps],
                        T( w[j][t], phases ) .== T( reg.vreg_pu * cj(reg.vreg_pu), phases )
                    )
                else  # use turn ratio
                    m[:kvl][j] = @constraint(m, [t in 1:net.Ntimesteps],
                        T( w[j][t], phases ) .== T( w[i][t], phases ) * reg.turn_ratio^2 
                    )
                end
                
            end
        end
    end
    nothing
end
