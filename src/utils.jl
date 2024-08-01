
function get_constraints_by_variable_name(m::JuMP.AbstractModel, v::AbstractString)
    ac = ConstraintRef[]
    for tup in list_of_constraint_types(m)
        append!(ac, all_constraints(m, tup[1], tup[2]))
    end
    filter( cr -> occursin(v, string(cr)), ac )
end


"""
    cj(A)

short cut for conj(transpose(A))
"""
function cj(A)
    conj(transpose(A))
end


"""
    phi_ij(j::String, net::Network, M::AbstractMatrix)

Down-select the matrix M by the phase from i -> j
"""
function phi_ij(j::String, net::Network, M::AbstractMatrix)
    N = convert(Matrix{GenericAffExpr{ComplexF64, VariableRef}}, [0 0im 0im; 0im 0. 0im; 0im 0im 0])
    for x in phases_into_bus(net, j), y in phases_into_bus(net, j)
        N[x,y] = M[x,y]
    end
    return N
end


"""
    phi_ij(j::String, net::Network, v::AbstractVector)

Down-select the vector v by the phase from i -> j
"""
function phi_ij(j::String, net::Network, v::AbstractVector)
    n = convert(Vector{GenericAffExpr{ComplexF64, VariableRef}}, [0im; 0im; 0im])
    for x in phases_into_bus(net, j)
        n[x] = v[x]
    end
    return n
end
