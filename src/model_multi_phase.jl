"""
    build_model!(m::JuMP.AbstractModel, p::Inputs{MultiPhase})

Add variables and constraints to `m` using the values in `p`. Calls the following functions:
```julia
add_variables(m, p)
constrain_power_balance(m, p)
constrain_substation_voltage(m, p)
constrain_KVL(m, p)
constrain_loads(m, p)
if p.relaxed
    constrain_cone(m, p)
else
    constrain_psd(m, p)
end
```
"""
function build_model!(m::JuMP.AbstractModel, p::Inputs{MultiPhase})

    add_variables(m, p)
    constrain_power_balance(m, p)

    constrain_KVL(m, p)
    constrain_loads(m, p)
    if p.relaxed
        constrain_cone(m, p)
    else
        constrain_psd(m, p)
    end
end

# hack for zero problem in MutableArithmetics
import Base: zero
zero(::Type{Any}) = 0.0
zero(::Type{Union{Float64, GenericAffExpr}}) = 0.0

# NEXT need to put zeros in matrices s.t. all variables are 3x1 or 3x3?
function add_variables(m, p::Inputs{MultiPhase})
    T = 1:p.Ntimesteps

    # map for phases into bus
    d = deepcopy(p.phases_into_bus)
    d[p.substation_bus] = [1,2,3]

    S = Dict{String, Matrix{Union{Float64, GenericAffExpr}}}
    # voltage squared is Hermitian, we make a dict with time and bus keys
    w = Dict{Int64, Dict{String, Matrix{Union{Float64, GenericAffExpr, ComplexF64}}}}()
    # current squared is Hermitian, we make a dict with time and bus keys
    l = Dict{Int64, S}()
    # complex line powers (at the sending end)
    Sij = Dict{Int64, S}()
    # complex powers injections
    Sj = Dict{Int64, S}()

    for t in T
        # fix head voltage at p.v0
        w[t] = Dict(p.substation_bus => [1.0+0im 0 0; 0 1 0; 0 0 1])
        # empty dicts for line values in each time step to fill
        l[t] = Dict{String, S}()
        Sij[t] = Dict{String, S}()
        # slack bus power injection
        Sj[t] = Dict(p.substation_bus =>  @variable(m, [1:3, 1:3] in ComplexPlane(), 
            base_name="Sj_" * string(t) *"_"* p.substation_bus,
            lower_bound = p.P_lo_bound + p.Q_lo_bound*im, 
            upper_bound = p.P_up_bound + p.Q_up_bound*im
        ))
        
        for j in keys(d)
            if j == p.substation_bus
                # already filled Sj and w and there are no lines into substation_bus
                continue
            end

            i = i_to_j(j, p)[1]
            i_j = string(i * "-" * j)  # for radial network there is only one i in i_to_j

            Sij[t][i_j] = convert(Matrix{Union{GenericAffExpr, Float64}}, zeros(3,3))
            l[t][i_j] = convert(Matrix{Union{GenericAffExpr, Float64}}, zeros(3,3))
            Sj[t][j] = convert(Matrix{Union{GenericAffExpr, Float64}}, zeros(3,1))
            w[t][j] = convert(Matrix{Union{GenericAffExpr, Float64}}, zeros(3,1))

            # fill in variables for all combinations of phases
            # TODO better bounds
            for phs1 in d[j], phs2 in d[j]
                Sij[t][i_j][phs1, phs2] = @variable(m, 
                    set = ComplexPlane(), base_name="Sij_" * string(t) *"_"* i_j *"_"*  string(phs1) * string(phs2), 
                    lower_bound = p.P_lo_bound + p.Q_lo_bound*im,
                    upper_bound = p.P_up_bound + p.Q_up_bound*im
                )
                l_up_bound = p.Isqaured_up_bounds[get_ijlinecode(i,j,p)]
                l[t][i_j][phs1, phs2] = @variable(m, 
                    set = ComplexPlane(), base_name="l_" * string(t) *"_"* i_j *"_"*  string(phs1) * string(phs2), 
                    lower_bound = 0.0 + 0.0*im,
                    upper_bound = l_up_bound + l_up_bound*im
                )
            end
            for phs in d[j]
                Sj[t][j][phs] = @variable(m, 
                    set = ComplexPlane(), base_name="Sj_" * string(t) *"_"* j *"_"* string(phs), 
                    lower_bound = p.P_lo_bound + p.Q_lo_bound*im,
                    upper_bound = p.P_up_bound + p.Q_up_bound*im
                )
                w[t][j][phs] = @variable(m, 
                    set = ComplexPlane(), base_name="w_" * string(t) *"_"* j *"_"* string(phs), 
                    lower_bound = p.v_lolim^2 + p.v_lolim^2*im,
                    upper_bound = p.v_uplim^2 + p.v_uplim^2*im
                )
            end
            # TODO line limit inputs
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
                Sⱼ[t][j] .- sum( diag(Sᵢⱼ[t][string(j*"-"*k)]) for k in j_to_k(j, p) ) .== 0
            )
        elseif isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # unconnected nodes
            @warn "Bus $j has no edges, setting Sⱼ and Qⱼ to zero."
            con = @constraint(m, [t in 1:p.Ntimesteps],
                diag(Sⱼ[t][j]) .== 0
            )
        elseif !isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # leaf nodes / sinks, flows in = draw out
            con = @constraint(m, [t in 1:p.Ntimesteps],
                sum( diag(Sᵢⱼ[t][string(i*"-"*j)]) for i in i_to_j(j, p) ) .-
                sum( diag(zij(i,j,p) * lᵢⱼ[t][string(i*"-"*j)]) for i in i_to_j(j, p) ) 
                + Sⱼ[t][j] .== 0
            )
        else
            con =  @constraint(m, [t in 1:p.Ntimesteps],
                sum( diag(Sᵢⱼ[t][string(i*"-"*j)]) for i in i_to_j(j, p) ) .- 
                sum( diag(zij(i,j,p) * lᵢⱼ[t][string(i*"-"*j)]) for i in i_to_j(j, p) ) .+
                Sⱼ[t][j] .- 
                sum( diag(Sᵢⱼ[t][string(j*"-"*k)]) for k in j_to_k(j, p) ) .== 0
            )
        end
        m[:loadbalcons][j] = con
    end
    # TODO Farivar and Low have b*v and q*v in these equations for shunts? neglected for now

    nothing
end


function cj(A)
    conj(transpose(A))
end


function constrain_KVL(m, p::Inputs{MultiPhase})
    w = m[:w]
    Sᵢⱼ = m[:Sij]
    lᵢⱼ = m[:l]
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            i_j = string(i*"-"*j)
            z = zij(i,j,p)
            # TODO only need the upper triangle of these element-wise, matrix constraints
            @constraint(m, [t in 1:p.Ntimesteps],
                w[t][j] .== w[t][i]
                    - (Sᵢⱼ[t][i_j] * cj(z) + z * cj(Sᵢⱼ[t][i_j]))
                    + z * lᵢⱼ[t][i_j] * cj(z)
            );
        end
    end
    nothing
end


function constrain_cone(m, p::Inputs{MultiPhase})
    w = m[:vsqrd]
    P = m[:Pᵢⱼ]
    Q = m[:Qᵢⱼ]
    l = m[:lᵢⱼ]
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            i_j = string(i*"-"*j)
            # # equivalent but maybe not handled as well in JuMP ?
            # @constraint(m, [t in 1:p.Ntimesteps],
            #     w[i,t] * l[i_j, t] ≥ P[i_j,t]^2 + Q[i_j,t]^2
            # )
            @constraint(m, [t in 1:p.Ntimesteps], 
                [w[i,t]/2, l[i_j, t], P[i_j,t], Q[i_j,t]] in JuMP.RotatedSecondOrderCone()
            )
        end
    end
end


function constrain_psd(m, p::Inputs{MultiPhase})
    w = m[:vsqrd]
    P = m[:Pᵢⱼ]
    Q = m[:Qᵢⱼ]
    l = m[:lᵢⱼ]
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            i_j = string(i*"-"*j)
            # need to use PSD?
            @constraint(m, [t in 1:p.Ntimesteps],
                w[i,t] * l[i_j, t] == P[i_j,t]^2 + Q[i_j,t]^2
            )
        end
    end
end


"""
    constrain_loads(m, p::Inputs{MultiPhase})

- set loads to negative of Inputs.Pload, which are normalized by Sbase when creating Inputs
- keys of Pload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus
"""
function constrain_loads(m, p::Inputs{MultiPhase})
    Pⱼ = m[:Pⱼ]
    Qⱼ = m[:Qⱼ]
    
    for j in p.busses
        if j in keys(p.Pload)
            @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,t] == -p.Pload[j][t] / p.Sbase
            )
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,t] == 0
            )
        end
        if j in keys(p.Qload)
            @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == -p.Qload[j][t] / p.Sbase
            )
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == 0
            )
        end
    end
    nothing
end
