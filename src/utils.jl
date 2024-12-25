
function get_constraints_by_variable_name(m::JuMP.AbstractModel, v::AbstractString)
    ac = ConstraintRef[]
    for tup in list_of_constraint_types(m)
        append!(ac, all_constraints(m, tup[1], tup[2]))
    end
    filter( cr -> occursin(v, string(cr)), ac )
end
