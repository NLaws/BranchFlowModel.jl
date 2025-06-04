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
        "Models" => ["single_phase_models.md", "multi_phase_models.md"],
        "Variables" => "variables.md",
        "Constraints" => "constraints.md",
        "Methods" => "methods.md",
        "Math" => "math.md"
    ],
    warnonly = true,  # TODO rm this and fix all the docs
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/NLaws/BranchFlowModel.jl.git"
)
