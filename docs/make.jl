using Documenter, Mooring

user_pages = [
  "UserGuide.md",
  "Installation"=>"userguide/installation.md",
  "Points"=>"userguide/points.md",
  "Segments"=>"userguide/segments.md",
  "Lines"=>"userguide/lines.md",
  "Materials"=>"userguide/materials.md",
]

developer_pages = []

entities_pages = [
  "MooringPoint"=>"API/Entities/MooringPoint.md",
  "MooringSegment"=>"API/Entities/MooringSegment.md",
  "MooringLine"=>"API/Entities/MooringLine.md",
]

geometry_pages = [
  "MooringDiscreteModel"=>"API/Geometry/MooringDiscreteModel.md",
]

io_pages = [
  "ParameterHandler"=>"API/IO/ParameterHandler.md",
]

physics_pages = [
  "Materials"=>"API/Physics/Materials.md",
  "EnvironmentalConditions"=>"API/Physics/EnvironmentalConditions.md",
  "SeaBed"=>"API/Physics/SeaBed.md",
  "Drag"=>"API/Physics/Drag.md",
  "PointMotion"=>"API/Physics/PointMotion.md",
  "TangentialDiffCalculus"=>"API/Physics/TangentialDiffCalculus.md",
]

api_pages = [ 
  "Overview"=>"API/API.md",
  "Entities" => entities_pages,
  "Geometry" => geometry_pages,
  "Input/Output" => io_pages,
  "Physics" => physics_pages
]

pages = [
    "Home" => "index.md",
    "Getting Started" => "GettingStarted.md",
    "User Guide" => user_pages,
    "Developer Guide" => "DeveloperGuide.md",
    "API" => api_pages,
]

makedocs(
    sitename = "Mooring.jl",
    format = Documenter.HTML(
    size_threshold=nothing,
    prettyurls=true,
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