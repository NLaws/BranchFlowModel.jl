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


# TODO test this function (esp. for matching openDSS results for IEEE13 system)
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


"""
    info_max_rpu_xpu(p::Inputs)

report the maximum per-unit resistance and reactance values of the lines.
It is important that the rpu ans xpu values be << 1. See Chiang and Baran 2013:
A load flow solution with feasible voltage magnitude always exists and is unique when
1. V0 ≈ 1
2. loss values < 1
3. rpu, xpu << 1
"""
function info_max_rpu_xpu(p::Inputs)
    Rmax = maximum([rij(i,j,p) for (i,j) in p.edges])
    Xmax = maximum([xij(i,j,p) for (i,j) in p.edges])
    @info("Max. Rpu: $Rmax   Max Xpu: $Xmax")
    return Rmax, Xmax
end
