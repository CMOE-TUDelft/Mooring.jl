using Documenter, Mooring

user_pages = []

developer_pages = []

entities_pages = [
  "API/Entities/MooringPoint.md",
  "API/Entities/MooringSegment.md",
  "API/Entities/MooringLine.md",
]

geometry_pages = [
  "API/Geometry/MooringDiscreteModel.md",
]

physics_pages = [
  "API/Physics/EnvironmentalConditions.md",
  "API/Physics/SeaBed.md",
  "API/Physics/Drag.md",
  "API/Physics/PointMotion.md",
  "API/Physics/TangentialDiffCalculus.md",
]

api_pages = [ 
  "API/API.md", 
  "Input/Output"=>"API/IO/ParameterHandler.md",
  "Entities" => entities_pages,
  "Geometry" => geometry_pages,
  "Physics" => physics_pages
]

pages = [
    "Home" => "index.md",
    "Getting Started" => "GettingStarted.md",
    "User Guide" => "UserGuide.md",
    "Developer Guide" => "DeveloperGuide.md",
    "API" => "API/API.md",
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