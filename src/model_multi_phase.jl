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
    constrain_substation_voltage(m, p)
    constrain_KVL(m, p)
    constrain_loads(m, p)
    if p.relaxed
        constrain_cone(m, p)
    else
        constrain_psd(m, p)
    end
end


function add_variables(m, p::Inputs{MultiPhase})
    T = 1:p.Ntimesteps

    # map for phases into bus
    d = deepcopy(p.phases_into_bus)
    d[p.substation_bus] = [1,2,3]
    
    # bus injections: indexed by time int, bus string, phase int
    # NOTE the index order has changed from SinglePhase to align with how we must index the Hermitian matrices
    @variable(m, Sⱼ[T, b in keys(d), d[b]] in ComplexPlane(), 
        lower_bound = p.P_lo_bound + p.Q_lo_bound*im,
        upper_bound = p.P_up_bound + p.Q_up_bound*im,
    )
    # NOTE that Sⱼ is not a dict, so it indexes like m[:Sⱼ][1, "675", 3]

    # voltage squared is Hermitian, we make a dict with time and bus keys
    w = Dict{Int64, Dict{String, Any}}()
    # current squared is Hermitian, we make a dict with time and bus keys
    l = Dict{Int64, Dict{String, Any}}()
    # complex line powers (at the sending end)
    Sij = Dict{Int64, Dict{String, Any}}()
    for t in T
        # fix head voltage at p.v0
        w[t] = Dict(p.substation_bus => repeat([p.v0], length(d[p.substation_bus])))
        l[t] = Dict{String, Any}()
        Sij[t] = Dict{String, Any}()
        for j in keys(d)
            if j == p.substation_bus
                continue
            end
            nphases = length(d[j])
            i = i_to_j(j, p)[1]
            i_j = string(i * "-" * j)  # for radial network there is only one i in i_to_j
            if nphases == 1
                w[t][j] = reshape(
                    [@variable(m, base_name="w_"*j)],
                    1,1  # make a 1x1 Matrix
                )
                l[t][i_j] = reshape(
                    [@variable(m, base_name="l_"*i_j)],
                    1,1  # make a 1x1 Matrix
                )
                Sij[t][i_j] = reshape(
                    [@variable(m, 
                        set = ComplexPlane(), base_name="Sij_"*i_j, 
                        lower_bound = p.P_lo_bound + p.Q_lo_bound*im,
                        upper_bound = p.P_up_bound + p.Q_up_bound*im
                    )],
                    1,1  # make a 1x1 Matrix
                )
            else
                w[t][j] = @variable(m, [1:nphases, 1:nphases] in HermitianPSDCone(), base_name="w_"*j)
                l[t][i_j] = @variable(m, [1:nphases, 1:nphases] in HermitianPSDCone(), base_name="l_"*i_j)
                Sij[t][i_j] = @variable(m, [1:nphases, 1:nphases] in ComplexPlane(), base_name="Sij_"*i_j,
                    lower_bound = p.P_lo_bound + p.Q_lo_bound*im, 
                    upper_bound = p.P_up_bound + p.Q_up_bound*im
                )
            end
            # bounds on diagonals, i.e. voltage squared on each phase
            @constraint(m, 
                p.v_lolim^2 .<= real.(diag(w[t][j])) .<= p.v_uplim^2
            )
            @constraint(m, 
                0.0 .<= real.(diag(l[t][i_j])) .<= p.Isqaured_up_bounds[get_ijlinecode(i,j,p)]
            )
            # @constraint(m, 
            #     p.P_lo_bound + p.Q_lo_bound*im .<= Sij[t][i_j] .<= p.P_up_bound + p.Q_up_bound*im
            # )
            # TODO line limit inputs
        end
    end
    m[:w] = w
    m[:l] = l
    m[:Sij] = Sij
    
    nothing
end

