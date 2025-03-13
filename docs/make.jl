using Documenter, Mooring

pages = [
    "Home" => "index.md",
    "Environmental Conditions" => "EnvironmentalConditions.md",
    "Sea Bed" => "SeaBed.md",
    "Drag" => "Drag.md",
    "Point Motion" => "PointMotion.md",
    "Tangential Differential Calculus" => "TangentialDiffCalculus.md",
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