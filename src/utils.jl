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


"""
    rij(i::AbstractString, j::AbstractString, p::Inputs{SinglePhase})

The per-unit resistance of line i->j
"""
function rij(i::AbstractString, j::AbstractString, p::Inputs{SinglePhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    rmatrix = p.Zdict[linecode]["rmatrix"] * linelength / p.Zbase
    return rmatrix[1]
end


"""
    xij(i::AbstractString, j::AbstractString, p::Inputs{SinglePhase})

The per-unit reacttance of line i->j
"""
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


function get_constraints_by_variable_name(m::JuMP.AbstractModel, v::AbstractString)
    ac = ConstraintRef[]
    for tup in list_of_constraint_types(m)
        append!(ac, all_constraints(m, tup[1], tup[2]))
    end
    filter( cr -> occursin(v, string(cr)), ac )
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
    phi_ij(j::String, p::Inputs, M::AbstractMatrix)

Down-select the matrix M by the phase from i -> j
"""
function phi_ij(j::String, p::Inputs, M::AbstractMatrix)
    N = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
    for x in p.phases_into_bus[j], y in p.phases_into_bus[j]
        N[x,y] = M[x,y]
    end
    return N
end


"""
    check_unique_solution_conditions(p::Inputs)


report the maximum per-unit immpedance and load values.
See Chiang and Baran 2013:
A load flow solution with feasible voltage magnitude always exists and is unique when
1. V0 ≈ 1
2. loss values < 1
3. rpu, xpu << 1
"""
function check_unique_solution_conditions(p::Inputs)
    Rmax, Xmax = info_max_rpu_xpu(p)
    maxP, maxQ = info_max_Ppu_Qpu(p)

    if p.v0 > 1.09 || p.v0 < 0.91
        @warn "The substation voltage should be set close to one."
    end

    if maxP > 0.1 || maxQ > 0.1 
        @warn "The maximum load is greater than 0.1 * Sbase. You should probably increase Sbase."
    end
    if Rmax > 1e-2 || Xmax > 1e-2
        @warn "\nThe per unit impedance values should be much less than one but the max is >0.01 .\
            \nYou should increase Zbase by increasing Vbase."
    end
    nothing
end


function reg_busses(p::Inputs)
    vcat( 
        getindex.(keys(p.regulators), 1),
        getindex.(keys(p.regulators), 2)
    )
end