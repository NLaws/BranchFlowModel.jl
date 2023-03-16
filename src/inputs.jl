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
        Isqaured_up_bounds::Dict{String, <:Real}
        phases_into_bus::Dict{String, Vector{Int}}
    end

# Inputs
- `edges` e.g. [("0", "1"), ("1", "2")]
- `linecodes` vector of string keys for the Zdict (impedance values for lines). When using an OpenDSS model a `linecode` is the `name` in `New linecode.name`
- `linelengths` vector of floats to scale impedance values
- `busses` vector of bus names
- `phases` vector of vectors for the line phases (for now just [[1], [1], ...])
- `Pload` dict with `busses` for keys and uncontrolled real power loads (positive is load)
- `Qload` dict with `busses` for keys and uncontrolled reactive power loads (positive is load)
- `Sbase` base apparent power for network, typ. feeder capacity. Used to normalize all powers
- `Vbase` base voltage for network, used to determine `Zbase`` = `Vbase`^2 / Sbase
- `Ibase` = `Sbase` / (`Vbase` * sqrt(3))
- `Zdict` dict with `linecodes` for keys and subdicts with "xmatrix" and "zmatrix" keys. Values are scaled by `Zbase`.
- `v0` slack bus reference voltage

TODO Zdict example
TODO test against simple model to make sure scaling is done right

!!! NOTE
    The `edges`, `linecodes`, `phases`, and `linelengths` are in mutual order (e.g. the i-th value in each list corresponds to the same line)
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
    Isqaured_up_bounds::Dict{String, <:Real}  # index on ij_edges = [string(i*"-"*j) for j in p.busses for i in i_to_j(j, p)]
    phases_into_bus::Dict{String, Vector{Int}}
    relaxed::Bool
    edge_keys::Vector{String}
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
        Isqaured_up_bounds=Dict{String, Float64}(),
        relaxed=true
    )

Lowest level Inputs constructor (the only one that returns the Inputs struct). 

!!! note
    The real and reactive loads provided are normalized using `Sbase`.
"""
function Inputs(
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
        Isqaured_up_bounds=Dict{String, Float64}(),
        relaxed=true
    )
    Ibase = Sbase / (Vbase * sqrt(3))
    # Ibase^2 should be used to recover amperage from lij ?
    Zbase = Vbase^2 / Sbase
    @info "Zbase: ", Zbase
    busses = String[]
    for t in edges
        push!(busses, t[1])
        push!(busses, t[2])
    end
    busses = unique(busses)

    if isempty(Isqaured_up_bounds)
        Isqaured_up_bounds = Dict(l => DEFAULT_AMP_LIMIT^2 for l in linecodes)
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
        Isqaured_up_bounds,
        phases_into_bus,
        relaxed,
        edge_keys
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
    d = open(dssfilepath) do io  # 
        parse_dss(io)  # method from PowerModelsDistribution
    end
    edges, linecodes, linelengths, linecodes_dict, phases, Isqaured_up_bounds = dss_dict_to_arrays(d, Sbase, Vbase)

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
        Isqaured_up_bounds=Isqaured_up_bounds,
        relaxed=relaxed
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