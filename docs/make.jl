using CommonOPF
using Documenter
using BranchFlowModel
using JuMP

makedocs(
    sitename = "BranchFlowModel",
    format = Documenter.HTML(),
    modules = [BranchFlowModel, CommonOPF],
    workdir = joinpath(@__DIR__, ".."),
    pages = [
        "User Documentation" => "index.md",
        "Methods" => "methods.md",
        "Decomposition" => "decomposition.md",
        "Math" => "math.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/NLaws/BranchFlowModel.jl.git"
)
