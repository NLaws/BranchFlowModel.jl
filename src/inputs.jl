"""
    mutable struct Inputs{T<:Phases} <: AbstractInputs
        edges::Array{Tuple, 1}
        linecodes::Array{String, 1}
        linelengths::Array{Float64, 1}
        busses::Array{String}
        phases::Vector{Vector}
        substation_bus::String
        Pload::Dict{String, Any}
        Qload::Dict{String, Any}
        Sbase::Real
        Vbase::Real
        Ibase::Real
        Zdict::Dict{String, Dict{String, Any}}
        v0::Real
        v_lolim::Real
        v_uplim::Real
        Zbase::Real
        Ntimesteps::Int
        pf::Float64
        Nnodes::Int
        P_up_bound::Float64
        Q_up_bound::Float64
        P_lo_bound::Float64
        Q_lo_bound::Float64
        Isquared_up_bounds::Dict{String, <:Real}
        phases_into_bus::Dict{String, Vector{Int}}
    end

Inputs
- `edges` Vector{Tuple} e.g. `[("0", "1"), ("1", "2")]`
- `linecodes` vector of string keys for the Zdict (impedance values for lines). When using an OpenDSS model a `linecode` is the `name` in `New linecode.name`
- `linelengths` vector of floats to scale impedance values
- `busses` vector of bus names
- `phases` vector of vectors with ints for the line phases (e.g. `[[1,2,3], [1,3], ...]`)
- `Pload` dict with `busses` for keys and uncontrolled real power loads (positive is load) by phase and time
- `Qload` dict with `busses` for keys and uncontrolled reactive power loads (positive is load) by phase and time
- `Sbase` base apparent power for network, typ. feeder capacity. Used to normalize all powers in model
- `Vbase` base voltage for network, used to determine `Zbase = Vbase^2 / Sbase`
- `Ibase` = `Sbase / (Vbase * sqrt(3))`
- `Zdict` dict with `linecodes` for keys and subdicts with "xmatrix" and "zmatrix" keys with per unit length values. Values are divided by `Zbase` and multiplied by linelength in mathematical model.
- `v0` slack bus reference voltage

TODO Zdict example

TODO test against simple model to make sure scaling is done right

!!! note
    The `edges`, `linecodes`, `phases`, `edge_keys`, and `linelengths` are in mutual order (i.e. the i-th value in each list corresponds to the same line)

"""
mutable struct Inputs{T<:Phases} <: AbstractInputs
    edges::Array{Tuple, 1}
    linecodes::Array{String, 1}
    linelengths::Array{Float64, 1}
    busses::Array{String}
    phases::Vector{Vector}
    substation_bus::String
    Pload::Dict{String, Any}
    Qload::Dict{String, Any}
    Sbase::Real
    Vbase::Real
    Ibase::Real
    Zdict::Dict{String, Dict{String, Any}}
    v0::Union{Real, AbstractVecOrMat{<:Number}}  # TODO MultiPhase v0 
    v_lolim::Real
    v_uplim::Real
    Zbase::Real
    Ntimesteps::Int
    pf::Float64
    Nnodes::Int
    P_up_bound::Float64
    Q_up_bound::Float64
    P_lo_bound::Float64
    Q_lo_bound::Float64
    Isquared_up_bounds::Dict{String, <:Real}  # index on ij_edges = [string(i*"-"*j) for j in p.busses for i in i_to_j(j, p)]
    phases_into_bus::Dict{String, Vector{Int}}
    relaxed::Bool
    edge_keys::Vector{String}
    regulators::Dict
end
# TODO line flow limits


"""
    Inputs(
        edges::Array{Tuple}, 
        linecodes::Array{String}, 
        linelengths::Array{Float64}, 
        phases::Vector{Vector},
        substation_bus::String;
        Pload, 
        Qload, 
        Sbase=1, 
        Vbase=1, 
        Zdict, 
        v0, 
        v_lolim=0.95, 
        v_uplim=1.05,
        Ntimesteps=1, 
        P_up_bound=1e4,
        Q_up_bound=1e4,
        P_lo_bound=-1e4,
        Q_lo_bound=-1e4,
        Isquared_up_bounds=Dict{String, Float64}(),
        relaxed=true
    )

Lowest level Inputs constructor (the only one that returns the Inputs struct). 

!!! note
    The real and reactive loads provided are normalized using `Sbase`.
"""
function Inputs(
        edges::AbstractVector{<:Tuple}, 
        linecodes::AbstractVector{<:AbstractString}, 
        linelengths::AbstractVector{<:Real}, 
        phases::AbstractVector{<:AbstractVector},
        substation_bus::String;
        Pload, 
        Qload, 
        Sbase=1, 
        Vbase=1, 
        Zdict, 
        v0, 
        v_lolim=0.95, 
        v_uplim=1.05,
        Ntimesteps=1, 
        P_up_bound=1e4,
        Q_up_bound=1e4,
        P_lo_bound=-1e4,
        Q_lo_bound=-1e4,
        Isquared_up_bounds=Dict{String, Float64}(),
        relaxed=true,
        regulators=Dict()
    )

    Ibase = Sbase / (Vbase * sqrt(3))
    # Ibase^2 should be used to recover amperage from lij ?
    Zbase = Vbase^2 / Sbase

    busses = String[]
    for t in edges
        push!(busses, t[1])
        push!(busses, t[2])
    end
    busses = unique(busses)

    if isempty(Isquared_up_bounds)
        Isquared_up_bounds = Dict(l => DEFAULT_AMP_LIMIT^2 for l in linecodes)
    end

    if v_lolim < 0 @error("lower voltage limit v_lolim cannot be less than zero") end
    if v_uplim < 0 @error("upper voltage limit v_uplim cannot be less than zero") end

    receiving_busses = collect(e[2] for e in edges)
    phases_into_bus = Dict(k=>v for (k,v) in zip(receiving_busses, phases))

    input_type = SinglePhase
    if any(get(v, "nphases", 1) > 1 for v in values(Zdict))
        input_type = MultiPhase
    end

    edge_keys = [string(i*"-"*j) for (i,j) in edges]

    Inputs{input_type}(
        edges,
        linecodes,
        linelengths,
        busses,
        phases,
        substation_bus,
        Pload,
        Qload,
        Sbase,
        Vbase,
        Ibase,
        Zdict,
        v0,
        v_lolim, 
        v_uplim,
        Zbase,
        Ntimesteps,
        0.1,  # power factor
        length(busses),  # Nnodes
        P_up_bound,
        Q_up_bound,
        P_lo_bound,
        Q_lo_bound,
        Isquared_up_bounds,
        phases_into_bus,
        relaxed,
        edge_keys,
        regulators
    )
end


"""
    Inputs(
        dssfilepath::String, 
        substation_bus::String;
        Pload::AbstractDict=Dict(), 
        Qload::AbstractDict=Dict(), 
        Sbase=1, 
        Vbase=1, 
        v0, 
        v_lolim=0.95, 
        v_uplim=1.05,
        Ntimesteps=1, 
        P_up_bound=1e4,
        Q_up_bound=1e4,
        P_lo_bound=-1e4,
        Q_lo_bound=-1e4,
        relaxed=true,
    )

Inputs constructor that parses a openDSS file for the network. If `Pload` and `Qload` are not provided
then the loads are also parsed from the openDSS file.
"""
function Inputs(
        dssfilepath::String, 
        substation_bus::String;
        Pload::AbstractDict=Dict(), 
        Qload::AbstractDict=Dict(), 
        Sbase=1.0, 
        Vbase=1.0, 
        v0=1.0, 
        v_lolim=0.95, 
        v_uplim=1.05,
        Ntimesteps=1, 
        P_up_bound=1e4,
        Q_up_bound=1e4,
        P_lo_bound=-1e4,
        Q_lo_bound=-1e4,
        relaxed=true,
    )

    d = let d 
        with_logger(SimpleLogger(Error)) do
            open(dssfilepath) do io  # 
                parse_dss(io)# method from PowerModelsDistribution
            end
        end
    end

    edges, linecodes, linelengths, linecodes_dict, phases, Isquared_up_bounds, regulators = dss_dict_to_arrays(d, Sbase, Vbase)

    if isempty(Pload) && isempty(Qload)
        Pload, Qload = dss_loads(d)
        # hack for single phase models
        if all(v == [1] for v in phases)
            # strip phase index out of loads
            newP = Dict{String, Any}()
            for (b,v) in Pload
                newP[b] = v[1]
            end
            Pload = newP
            newQ = Dict{String, Any}()
            for (b,v) in Qload
                newQ[b] = v[1]
            end
            Qload = newQ
        end
    end

    # TODO line limits from OpenDSS ?

    Inputs(
        edges,
        linecodes,
        linelengths,
        phases,
        substation_bus;
        Pload=Pload, 
        Qload=Qload,
        Sbase=Sbase, 
        Vbase=Vbase, 
        Zdict=linecodes_dict, 
        v0=v0,
        v_lolim = v_lolim, 
        v_uplim = v_uplim, 
        Ntimesteps=Ntimesteps,
        P_up_bound=P_up_bound,
        Q_up_bound=Q_up_bound,
        P_lo_bound=P_lo_bound,
        Q_lo_bound=Q_lo_bound,
        Isquared_up_bounds=Isquared_up_bounds,
        relaxed=relaxed,
        regulators=regulators
    )
end


"""
    singlephase38linesInputs(;
        Pload=Dict{String, AbstractArray{Real, 1}}(), 
        Qload=Dict{String, AbstractArray{Real, 1}}(), 
        T=24,
        loadnodes = ["3", "5", "36", "9", "10", "11", "12", "13", "15", "17", "18", "19", "22", "25", 
                    "27", "28", "30", "31", "32", "33", "34", "35"],
        Sbase = 1e6,
        Vbase = 12.5e3,
        v0=1.0,
        v_uplim = 1.05,
        v_lolim = 0.95,
    )

Convenience function for creating a single phase network with 38 lines and nodes. 
Taken from:
Andrianesis et al. 2019 "Locational Marginal Value of Distributed Energy Resources as Non-Wires Alternatives"

NOTE that Inputs is a mutable struct (s.t. loads can be added later).
"""
function singlephase38linesInputs(;
        Pload=Dict{String, AbstractArray{Real, 1}}(), 
        Qload=Dict{String, AbstractArray{Real, 1}}(), 
        T=24,
        loadnodes = ["3", "5", "36", "9", "10", "11", "12", "13", "15", "17", "18", "19", "22", "25", 
                    "27", "28", "30", "31", "32", "33", "34", "35"],
        Sbase = 1e6,
        Vbase = 12.5e3,
        v0=1.0,
        v_uplim = 1.05,
        v_lolim = 0.95,
    )

    if isempty(Pload)  # fill in default loadnodes
        Pload = Dict(k => Real[] for k in loadnodes)
    end
    if isempty(Qload)  # fill in default loadnodes
        Qload = Dict(k => Real[] for k in loadnodes)
    end

    Inputs(
        joinpath(dirname(@__FILE__), "..", "test", "data", "singlephase38lines", "master.dss"), 
        "0";
        Pload=Pload, 
        Qload=Qload,
        Sbase=Sbase, 
        Vbase=Vbase, 
        v0 = v0,
        v_uplim = v_uplim,
        v_lolim = v_lolim,
        Ntimesteps = T
    )
end


"""
    delete_edge_ij!(i::String, j::String, p::Inputs{BranchFlowModel.SinglePhase})

delete edge `(i, j)` from
- p.edges
- p.phases
- p.linelengths
- p.edge_keys
- p.Isquared_up_bounds

NOTE do not delete!(p.Zdict, ij_linecode) nor delete!(p.Isquared_up_bounds, ij_linecode) 
because anything indexed on linecodes can be used for multiple lines
"""
function delete_edge_ij!(i::String, j::String, p::Inputs{BranchFlowModel.SinglePhase})
    idx = get_ij_idx(i, j, p)
    deleteat!(p.edges,       idx)
    deleteat!(p.linecodes,   idx)
    deleteat!(p.phases,      idx)
    deleteat!(p.linelengths, idx)
    deleteat!(p.edge_keys,   idx)
    true
end


"""
    delete_bus_j!(j::String, p::Inputs{BranchFlowModel.SinglePhase})

Remove bus `j` from `p.busses`
"""
function delete_bus_j!(j::String, p::Inputs{BranchFlowModel.SinglePhase})
    p.busses = setdiff(p.busses, [j])
    if j in keys(p.Pload)
        delete!(p.Pload, j)
    end
    if j in keys(p.Qload)
        delete!(p.Qload, j)
    end
    true
end


"""
    reduce_tree!(p::Inputs{SinglePhase})

combine any line sets with intermediate busses that have indegree == outdegree == 1
and is not a load bus into a single line
"""
function reduce_tree!(p::Inputs{BranchFlowModel.SinglePhase})
    # TODO make graph once in Inputs ?
    g = make_graph(p.busses, p.edges)
    int_bus_map = get_prop(g, :int_bus_map)
    reducable_buses = String[]
    load_buses = Set(vcat(collect(keys(p.Pload)), collect(keys(p.Qload))))
    for v in vertices(g)
        if indegree(g, v) == outdegree(g, v) == 1 && !(int_bus_map[v] in load_buses)
            push!(reducable_buses, int_bus_map[v])
        end
    end
    @debug("Removing the following busses: \n$reducable_buses")
    # replace two lines with one
    for j in reducable_buses
        # get all the old values
        i, k = i_to_j(j, p)[1], j_to_k(j, p)[1]
        ij_idx, jk_idx = get_ij_idx(i, j, p), get_ij_idx(j, k, p)
        ij_len, jk_len = p.linelengths[ij_idx], p.linelengths[jk_idx]
        ij_linecode, jk_linecode = get_ijlinecode(i,j,p), get_ijlinecode(j,k,p)
        r_ij, x_ij, r_jk, x_jk = rij(i,j,p)*p.Zbase, xij(i,j,p)*p.Zbase, rij(j,k,p)*p.Zbase, xij(j,k,p)*p.Zbase
        # make the new values
        r_ik = r_ij + r_jk
        x_ik = x_ij + x_jk
        ik_len = ij_len + jk_len
        ik_linecode = ik_key = i * "-" * k
        ik_amps = minimum([p.Isquared_up_bounds[ij_linecode], p.Isquared_up_bounds[jk_linecode]])
        # delete the old values
        delete_edge_ij!(i, j, p)
        delete_edge_ij!(j, k, p)
        delete_bus_j!(j, p)
        # add the new values
        push!(p.edges, (i, k))
        push!(p.linecodes, ik_linecode)
        push!(p.phases, [1])
        push!(p.linelengths, ik_len)
        push!(p.edge_keys, ik_key)
        p.Zdict[ik_linecode] = Dict(
            "nphases" => 1,
            "name" => ik_linecode,
            "rmatrix" => [r_ik / ik_len],
            "xmatrix" => [x_ik / ik_len],
        )
        p.Isquared_up_bounds[ik_linecode] = ik_amps
    end
    @info("Removed $(length(reducable_buses)) busses.")
end


function trim_tree_once!(p::Inputs{BranchFlowModel.SinglePhase})
    trimmable_busses = setdiff(leaf_busses(p), union(keys(p.Pload), keys(p.Qload)))
    if isempty(trimmable_busses) return false end
    trimmable_edges = Tuple[]
    for j in trimmable_busses
        for i in i_to_j(j, p)
            push!(trimmable_edges, (i,j))
        end
    end
    @debug("Deleting the following edges from the Inputs:")
    for edge in trimmable_edges @debug(edge) end
    for (i,j) in trimmable_edges
        delete_edge_ij!(i, j, p)
        delete_bus_j!(j, p)
    end
    true
end


function trim_tree!(p::Inputs{BranchFlowModel.SinglePhase})
    n_edges_before = length(p.edges)
    trimming = trim_tree_once!(p)
    while trimming
        trimming = trim_tree_once!(p)
    end
    n_edges_after = length(p.edges)
    @info("Removed $(n_edges_before - n_edges_after) edges.")
    true
end


function make_sub_inputs(
    p::Inputs{BranchFlowModel.SinglePhase}, 
    edges_to_delete::Vector, 
    busses_to_delete::Vector{String}
    )
    pc = deepcopy(p)
    for e in edges_to_delete
        delete_edge_ij!(e[1], e[2], pc)
    end
    for b in busses_to_delete
        delete_bus_j!(b, pc)
    end
    return pc
end


"""
    split_inputs(p::Inputs{BranchFlowModel.SinglePhase}, bus::String, g::SimpleDiGraph)

Split inputs into one graph for everything above `bus` and one graph for everything
    below `bus`.
"""
function split_inputs(p::Inputs{BranchFlowModel.SinglePhase}, bus::String)
    g = make_graph(p.busses, p.edges)
    in_buses = collect(all_inneighbors(g, bus, String[]))
    out_buses = collect(all_outneighbors(g, bus, String[], String[]))
    # in/out_buses do not have bus, but sub_busses does have bus
    # we want to keep bus in both Inputs

    sub_busses, sub_edges = induced_subgraph(g, vcat(out_buses, bus))
    p_above = make_sub_inputs(p, sub_edges, out_buses)

    sub_busses, sub_edges = induced_subgraph(g, vcat(in_buses, bus))
    p_below = make_sub_inputs(p, sub_edges, in_buses)
    p_below.substation_bus = bus

    return p_above, p_below
end



"""
    split_inputs(p::Inputs{BranchFlowModel.SinglePhase}, bus::String, out_buses::Vector{String})

Split `p` into `p_above` and `p_below` where `p_below` has only `out_buses` and `p_above` has
`union( [bus], setdiff(p.busses, out_buses) )`.

Note that `out_buses` must contain `bus`
"""
function split_inputs(p::Inputs{BranchFlowModel.SinglePhase}, bus::String, out_buses::Vector{String})
    if !(bus in out_buses)
        throw(@error "Cannot split inputs: bus is not in out_buses.")
    end
    g = make_graph(p.busses, p.edges)
    in_buses = setdiff(p.busses, out_buses)
    # in_buses does not have bus, but sub_busses does have bus
    # we want to keep bus in both p_above and p_below

    sub_busses, edges_to_remove = induced_subgraph(g, out_buses)
    # NOTE edges_to_remove does not include edges going out of out_buses
    p_above = make_sub_inputs(p, edges_to_remove, setdiff(out_buses, [bus]))

    sub_busses, edges_to_remove = induced_subgraph(g, vcat(in_buses, bus))
    p_below = make_sub_inputs(p, edges_to_remove, in_buses)
    p_below.substation_bus = bus

    return p_above, p_below
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


function info_max_Ppu_Qpu(p::Inputs)
    maxP = maximum(maximum.(values(p.Pload))) / p.Sbase
    maxQ = maximum(maximum.(values(p.Qload))) / p.Sbase
    @info("Max. Rpu: $maxP   Max Xpu: $maxQ")
    return maxP, maxQ
end