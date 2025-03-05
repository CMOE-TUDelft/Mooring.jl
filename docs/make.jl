using Documenter, Mooring

pages = [
    "Home" => "index.md",
    "Mooring.StressLinear" => "StressLinear.md"
]

makedocs(
    sitename = "Mooring.jl",
    remotes = nothing,
    pages = pages
)