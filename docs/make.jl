using Documenter, Mooring

pages = [
    "Home" => "index.md",
    "StressLinear" => "StressLinear.md"
]

makedocs(
    sitename = "Mooring.jl",
    format = Documenter.HTML(
    size_threshold=nothing
    ),
    modules = [Mooring],
    doctest = false,
    warnonly = [:cross_references,:missing_docs],
    checkdocs = :exports,
    remotes = nothing,
    pages = pages
)

deploydocs(
  repo = "github.com/CMOE-TUDelft/Mooring.jl.git",
)