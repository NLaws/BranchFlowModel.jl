"""
    function i_to_j(j::AbstractString, p::Inputs)
find all busses upstream of bus j

!!! note
    In a radial network this function should return an Array with length of 1.
"""
function i_to_j(j::AbstractString, p::Inputs)
    convert(Array{String, 1}, map(x->x[1], filter(t->t[2]==j, p.edges)))
end


"""
    function j_to_k(j::AbstractString, p::Inputs)
find all busses downstream of bus j
"""
function j_to_k(j::AbstractString, p::Inputs)
    convert(Array{String, 1}, map(x->x[2], filter(t->t[1]==j, p.edges)))
end


function rij(i::AbstractString, j::AbstractString, p::Inputs{SinglePhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    rmatrix = p.Zdict[linecode]["rmatrix"] * linelength / p.Zbase
    return rmatrix[1]
end


function xij(i::AbstractString, j::AbstractString, p::Inputs{SinglePhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    xmatrix = p.Zdict[linecode]["xmatrix"] * linelength / p.Zbase
    return xmatrix[1]
end


function rij(i::AbstractString, j::AbstractString, p::Inputs{MultiPhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    rmatrix = p.Zdict[linecode]["rmatrix"] * linelength / p.Zbase
    return rmatrix
end


function xij(i::AbstractString, j::AbstractString, p::Inputs{MultiPhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    xmatrix = p.Zdict[linecode]["xmatrix"] * linelength / p.Zbase
    return xmatrix
end


function zij(i::AbstractString, j::AbstractString, p::Inputs{MultiPhase})
    r = rij(i, j, p)
    x = xij(i, j, p)
    # need to fill out to 3x3 with zeros
    phases = sort(p.phases_into_bus[j])
    z = convert(Matrix{ComplexF64}, zeros(3,3))
    for (ii, phs1) in enumerate(phases), (jj, phs2) in enumerate(phases)
        z[phs1, phs2] = r[ii, jj] + x[ii, jj]*im
    end
    return z
end


function get_ijlinelength(i::AbstractString, j::AbstractString, p::Inputs)
    ij_idx = get_ij_idx(i, j, p)
    return p.linelengths[ij_idx]
end


function get_ijlinecode(i::AbstractString, j::AbstractString, p::Inputs)
    ij_idx = get_ij_idx(i, j, p)
    return p.linecodes[ij_idx]
end


function get_ijedge(i::AbstractString, j::AbstractString, p::Inputs)
    ij_idx = get_ij_idx(i, j, p)
    return p.edges[ij_idx]
end


function get_ij_idx(i::AbstractString, j::AbstractString, p::Inputs)
    ij_idxs = findall(t->(t[1]==i && t[2]==j), p.edges)
    if length(ij_idxs) > 1
        error("found more than one edge for i=$i and j=$j")
    elseif length(ij_idxs) == 0
        error("found no matching edges for i=$i and j=$j")
    else
        return ij_idxs[1]
    end
end


function get_edge_values(var_prefix::AbstractString, m::JuMP.AbstractModel, p::Inputs)
    vals = Float64[]
    for edge in p.edges
        var = string(var_prefix, "[", edge[1], "-", edge[2], "]")
        try
            val = value(variable_by_name(m, var))
            if startswith(var_prefix, "l")
                val = sqrt(val)
            end
            push!(vals, round(val; digits=8))
        catch e
            println(var, "failed", e)
        end
    end
    return vals
end


function get_bus_values(var::Symbol, m::JuMP.AbstractModel, p::Inputs{SinglePhase})
    vals = value.(m[var])
    d = Dict()
    for b in p.busses
        d[b] = vals[b,:].data
    end
    return d
end


function get_constraints_by_variable_name(m::JuMP.AbstractModel, v::AbstractString)
    ac = ConstraintRef[]
    for tup in list_of_constraint_types(m)
        append!(ac, all_constraints(m, tup[1], tup[2]))
    end
    filter( cr -> occursin(v, string(cr)), ac )
end


"""
    check_soc_inequalities(m::JuMP.AbstractModel, p::Inputs)

create and return a vector of the gaps in the second order cone constraints
"""
function check_soc_inequalities(m::JuMP.AbstractModel, p::Inputs)
    @assert termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]
    w = value.(m[:vsqrd])
    lij = value.(m[:lij])
    Pj = value.(m[:Pⱼ])
    Qj = value.(m[:Qⱼ])
    gaps = Vector{Float64}()  # todo size empty array based on number of SOC cons

    # TODO mesh accounting
    for j in p.busses, i in i_to_j(j, p), t in 1:p.Ntimesteps
        push!(gaps, 
            lij[string(i*"-"*j), t] * w[i, t] - Pj[j,t]^2 + Qj[j,t]^2 
        )
    end
    return gaps
end


"""
    get_load_bal_shadow_prices(m::JuMP.AbstractModel, p::Inputs)

create and return a dict indexed by bus and time for shadow prices
    (just real prices for now)
"""
function get_load_bal_shadow_prices(m::JuMP.AbstractModel, p::Inputs)
    @assert has_duals(m)
    d = Dict{String, Dict}(j => Dict{Int64, Float64}() for j in p.busses)

    for j in p.busses, t in 1:p.Ntimesteps
        d[j][t] = JuMP.shadow_price(m[:loadbalcons][j]["p"][t])
    end
    return d
end

"""
    cj(A)

short cut for conj(transpose(A))
"""
function cj(A)
    conj(transpose(A))
end


"""
    check_rank_one(m::JuMP.AbstractModel, p::Inputs{BranchFlowModel.MultiPhase}, tol=1e-3)

Check the rank of the `m[:H]` matrices from the PSD cone constraints.
Warnings express any values with rank greater than one.
"""
function check_rank_one(m::JuMP.AbstractModel, p::Inputs{BranchFlowModel.MultiPhase}, tol=1e-3)
    for j in p.busses, t in 1:p.Ntimesteps
        if j == p.substation_bus continue end
        eigs = eigvals!(JuMP.value.(m[:H][t][j]))
        eigs ./= maximum(eigs)
        rank_H = sum(map(x-> x > tol ? 1 : 0, eigs))
        if rank_H > 1
            @warn("Bus $j in time step $t has H matrix of rank $rank_H")
        end
    end
end


"""
    phi_ij(i::String, j::String, p::Inputs, M::AbstractMatrix)

Down-select the matrix M by the phase from i -> j
"""
function phi_ij(i::String, j::String, p::Inputs, M::AbstractMatrix)
    N = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
    for x in p.phases_into_bus[j], y in p.phases_into_bus[j]
        N[x,y] = M[x,y]
    end
    return N
end

