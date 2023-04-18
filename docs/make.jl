using Documenter
using BranchFlowModel
using JuMP

makedocs(
    sitename = "BranchFlowModel",
    format = Documenter.HTML(),
    modules = [BranchFlowModel],
    workdir = joinpath(@__DIR__, ".."),
    pages = [
        "User Documentation" => "index.md",
        "Methods" => "methods.md",
        "Math" => "math.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/NLaws/BranchFlowModel.jl.git"
)
