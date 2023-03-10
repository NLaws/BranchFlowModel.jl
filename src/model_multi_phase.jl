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
        v0 = [1.0+0im; -0.5 + im*sqrt(3)/2; -0.5 - im*sqrt(3)/2]
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

                        l[t][i_j][phs1, phs2] = @variable(m, 
                            base_name="l_" * string(t) *"_"* i_j *"_"*  string(phs1) * string(phs2), 
                            lower_bound = 0.0,
                            upper_bound = l_up_bound
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
                            lower_bound = 0.0 + 0.0*im,
                            upper_bound = l_up_bound + l_up_bound*im
                        )

                        w[t][j][phs1, phs2] = @variable(m, 
                            set = ComplexPlane(), base_name="w_" * string(t) *"_"* j *"_"* string(phs1) * string(phs2), 
                            lower_bound = p.v_lolim^2 + p.v_lolim^2*im,
                            upper_bound = p.v_uplim^2 + p.v_uplim^2*im,
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


function matrix_upper_triangle_to_vec(M::AbstractMatrix{T}, phases::AbstractVector{Int}) where T
    v = T[]
    for i in phases, j in phases 
        if j >= i
            push!(v, M[i,j])
        end
    end
    return v
end


function constrain_KVL(m, p::Inputs{MultiPhase})
    w = m[:w]
    Sᵢⱼ = m[:Sij]
    lᵢⱼ = m[:l]

    T = matrix_upper_triangle_to_vec
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

- set loads to negative of Inputs.Pload, which are normalized by Sbase when creating Inputs
- keys of Pload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus
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
