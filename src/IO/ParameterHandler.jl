module ParameterHandler

# using YAML
# using JSON


"""
 PointParameters struct

  This struct is used to define a point in the mooring system. It contains all relevant information
  about the point's position, motion, and mesh size.
"""
struct PointParameters
    id::Int
    tag::String
    coords::Vector{Float64}
    motion::String
    mesh_size::Float64
end

"""
    default_point(id::Int; tag="Point_1", coords=[0.0,0.0,0.0], motion=nothing, mesh_size=1.0) -> PointParameters

Return a default `Point` struct with basic values.

# Arguments
- `id`: Point identifier
- `tag`: Human-readable name
- `coords`: Coordinates in 2D or 3D
- `motion`: Optional motion definition
- `mesh_size`: Characteristic mesh size for meshing

# Example
```julia
p = default_point(1; coords=[0.0, 10.0, -5.0])
```
"""
function default_point(id::Int; coords=[0.0,0.0,0.0], motion="default", mesh_size=1.0)
  tag = "Point_$id"
  return PointParameters(id, tag, coords, motion, mesh_size)
end

"""
  SegmentParameters

This struct is used to define a segment in the mooring system. It contains all relevant information
about the segment's start and stop points, length, material properties, density, and cross-sectional area.
"""
struct SegmentParameters
  id::Int
  tag::String
  start_point::Int
  stop_point::Int
  length::Float64
  material::String
  density::Float64
  area::Float64
end

"""
default_segment(id::Int; tag="S", start=1, stop=2, length=10.0, material="steel", density=7850.0, area=0.01) -> SegmentParameters

Return a default Segment struct with nominal values.
"""
function default_segment(id::Int; start=1, stop=2, length=10.0, material="default", density=7850.0, area=0.01)
  tag = "Segment_$id"
  return SegmentParameters(id, tag, start, stop, length, material, density, area)
end

"""
LineParameters

  This struct is used to define a line in the mooring system. It contains the list of points and segments
  of the line.
"""
struct LineParameters
  id::Int
  tag::String
  points::Vector{Int}
  segments::Vector{Int}
end

"""
default_line(id::Int; points=[1,2], segments=[1]) -> LineParameters

Return a default Line consisting of a list of point IDs and segment IDs.
"""
function default_line(id::Int; points=[1,2], segments=[1])
  tag = "Line_$id"
  return LineParameters(id, tag, points, segments)
end

"""
DragParameters

  This struct is used to define the drag properties of a segment in the mooring system.
  The default type is "NoDrag", which does not require additional parameters. Possible drag
  types are:
  - "NoDrag": No drag forces
  - "Custom": Custom drag force from the following custom coefficients:
      Required input parameters:
      - `nd::Real`: Nominal diameter [m]
      - `AStr::Real`: Area (structural) [m2]
      Optional parameters (with default values):
      - `Cd_n::Real`: Drag coefficient normal
      - `Cd_t::Real`: Drag coefficient tangent
      - `Ca_n::Real`: Added mass coefficient
      - `Ca_t::Real`: Added mass coefficient
      - `ρw::Real = 1025`: Water density [Kg/m3]
      - `od::Real = nd`: Outer diameter [m]
      - `id::Real = 0.0`: Inner diameter [m]
      - `dd_n::Real = nd`: Drag diameter normal
      - `dd_t::Real = nd / π`: Drag diameter tangent
  - "ChainStudless": Chain studless drag properties
      Required input parameters:
      - `nd::Real`: Nominal diameter [m]
      - `AStr::Real`: Area (structural) [m2]
      Optional parameters (with default values):
      - `Cd_n::Real = 2.4`: Drag coefficient normal
      - `Cd_t::Real = 1.15`: Drag coefficient tangent
      - `Ca_n::Real = 1.0`: Added mass coefficient
      - `Ca_t::Real = 0.5`: Added mass coefficient
      - `ρw::Real = 1025`: Water density [Kg/m3]
      - `od::Real = 1.80 * nd`: Outer diameter [m]
      - `id::Real = 0.0`: Inner diameter [m]
      - `dd_n::Real = nd`: Drag diameter normal
      - `dd_t::Real = nd / π`: Drag diameter tangent
"""
struct DragParameters
  dragType::String
  ρw::Float64
  nd::Float64
  od::Float64
  id::Float64
  AStr::Float64
  Cd_n::Float64
  Cd_t::Float64
  dd_n::Float64
  dd_t::Float64
  Ca_n::Float64
  Ca_t::Float64
  Cfd_n::Float64
  Cfd_t::Float64
end
"""
default_drag(; dragType="NoDrag", ρw=0.0, nd=0.0, od=0.0, id=0.0, AStr=0.0,
Cd_n=0.0, Cd_t=0.0, dd_n=0.0, dd_t=0.0,
Ca_n=0.0, Ca_t=0.0, Cfd_n=0.0, Cfd_t=0.0) -> DragProperties

Return default parameters for NoDrag option.
"""
function default_drag(; dragType="NoDrag", ρw=0.0, nd=0.1, od=0.1, id=0.0, AStr=0.01,
  Cd_n=0.0, Cd_t=0.0, dd_n=0.0, dd_t=0.0,
  Ca_n=0.0, Ca_t=0.0, Cfd_n=0.0, Cfd_t=0.0)
  return DragParameters(dragType, ρw, nd, od, id, AStr, Cd_n, Cd_t, dd_n, dd_t, Ca_n, Ca_t, Cfd_n, Cfd_t)
end


#### HERE!

"""
default_waves(; Hs=0.0, Tp=0.0, h0=100.0, nω=64, seed=0, ωc=-1.0, enableWaveSpec=false) -> WaveParameters

Return default wave parameters (no waves by default).
"""
function default_waves(; Hs=0.0, Tp=0.0, h0=100.0, nω=64, seed=0, ωc=-1.0, enableWaveSpec=false)
  return WaveParameters(Hs, Tp, h0, nω, seed, ωc, enableWaveSpec)
end

"""
default_material(name::String="steel"; E=2.1e11, mu=0.3) -> Material

Return a default Material struct with elastic properties.
"""
function default_material(name::String="steel"; E=2.1e11, mu=0.3)
  return Material(name, Dict("E"=>E, "mu"=>mu))
end

"""
default_motion(id::Int=1; type="fixed", properties=Dict()) -> Motion

Return a default motion (fixed point by default).
"""
function default_motion(id::Int=1; type="fixed", properties=Dict())
  return Motion(type, properties)
end

"""
default_seabed(; kn=30e3, linear_damping_factor=0.05, quadratic_damping_factor=0.0,
od=0.1, A=0.008, tanh_ramp=1.0e2, penetration_depth_ramp=1.0e-3,
still_weight=0.0) -> SeaBedParameters

Return default seabed parameters with moderate stiffness and light damping.
"""
function default_seabed(; kn=30e3, linear_damping_factor=0.05, quadratic_damping_factor=0.0,
  od=0.1, A=0.008, tanh_ramp=1.0e2, penetration_depth_ramp=1.0e-3,
  still_weight=0.0)
  cnstz = kn * od / A
  return SeaBedParameters(kn, linear_damping_factor, quadratic_damping_factor,
  od, A, tanh_ramp, penetration_depth_ramp,
  still_weight, cnstz)
end

---------------------------------------------------------------------
Central ParameterHandler Struct
---------------------------------------------------------------------
"""
ParameterHandler

Unified container for all experiment parameters in the Mooring library.

Fields:

points::Dict{Int, Point}

segments::Dict{Int, Segment}

lines::Dict{Int, Line}

drag::Dict{String, DragProperties}

waves::WaveParameters

materials::Dict{String, Material}

motions::Dict{Int, Motion}

seabed::SeaBedParameters
"""
struct ParameterHandler
  points::Dict{Int, Point}
  segments::Dict{Int, Segment}
  lines::Dict{Int, Line}
  drag::Dict{String, DragProperties}
  waves::WaveParameters
  materials::Dict{String, Material}
  motions::Dict{Int, Motion}
  seabed::SeaBedParameters
end

"""
from_defaults() -> ParameterHandler

Construct a ParameterHandler with a single point, segment, line, drag property,
material, wave definition, and seabed parameters using the default_* functions.
"""
function from_defaults()
  return ParameterHandler(
  Dict(1 => default_point(1)),
  Dict(1 => default_segment(1)),
  Dict(1 => default_line(1)),
  Dict("default" => default_drag()),
  default_waves(),
  Dict("steel" => default_material("steel")),
  Dict(1 => default_motion(1)),
  default_seabed()
  )
end

"""
  ParameterHandler

A unified container for all experiment parameters in the Mooring library.
It stores the following fields:

# Fields
- `points::Dict{Int, Point}`: Dictionary of topological points by ID.
- `segments::Dict{Int, Segment}`: Dictionary of segments by ID.
- `lines::Dict{Int, Line}`: Dictionary of lines by ID.
- `drag::Dict{String, DragProperties}`: Dictionary of drag property sets.
- `waves::WaveParameters`: Global wave parameters.
- `materials::Dict{String, Material}`: Dictionary of materials.
- `motions::Dict{Int, Motion}`: Dictionary of motions for points.
- `seabed::SeaBedParameters`: Global seabed interaction parameters.

All input methods (`from_defaults`, `from_yaml`, `from_json`, `from_dat_folder`)
produce a consistent `ParameterHandler` object.
"""
struct ParameterHandler
  points::Dict{Int, Point}
  segments::Dict{Int, Segment}
  lines::Dict{Int, Line}
  drag::Dict{String, DragProperties}
  waves::WaveParameters
  materials::Dict{String, Material}
  motions::Dict{Int, Motion}
  seabed::SeaBedParameters
end


# """
#     from_defaults() -> ParameterHandler

# Initialize a ParameterHandler with default values.
# """
# function from_defaults()
#     return ParameterHandler(
#         Dict{Int, Point}(),
#         Dict{Int, Segment}(),
#         Dict{Int, Line}(),
#         Dict{String, DragProperties}(),
#         WaveParameters(),  # default constructor
#         Dict{String, Material}(),
#         Dict{Int, Motion}(),
#         SeaBedParameters()
#     )
# end

# """
#     from_yaml(path::String) -> ParameterHandler

# Parse YAML input into a ParameterHandler.
# """
# function from_yaml(path::String)
#     data = YAML.load_file(path)
#     return parse_dict(data)
# end

# """
#     from_json(path::String) -> ParameterHandler
# """
# function from_json(path::String)
#     data = JSON.parsefile(path)
#     return parse_dict(data)
# end

# """
#     from_dat_folder(folder::String) -> ParameterHandler

# Aggregate multiple .dat files into a single ParameterHandler and 
# produce a YAML export.
# """
# function from_dat_folder(folder::String)
#     # (placeholder) implement file aggregation
#     data = Dict()  # aggregated dictionary from .dat files
#     handler = parse_dict(data)

#     # export aggregated version for reproducibility
#     to_yaml(handler, joinpath(folder, "aggregated.yaml"))
#     return handler
# end

# # --- Export ---
# function to_yaml(ph::ParameterHandler, path::String)
#     open(path, "w") do io
#         YAML.write(io, to_dict(ph))
#     end
# end

# function to_json(ph::ParameterHandler, path::String)
#     open(path, "w") do io
#         JSON.print(io, to_dict(ph))
#     end
# end

# # --- Internal Conversion Utilities ---
# function parse_dict(data::Dict)::ParameterHandler
#     # Map Dict to strongly typed structs
#     points = Dict(id => Point(id, d["tag"], d["coordinates"], d["motion"], d["mesh_size"]) for (id, d) in data["points"])
#     segments = Dict(id => Segment(id, d["tag"], d["start"], d["stop"], d["length"],
#                                   d["material"], d["density"], d["area"]) for (id, d) in data["segments"])
#     lines = Dict(id => Line(d["points"], d["segments"]) for (id, d) in data["lines"])
#     drag = Dict(name => DragProperties(; d... ) for (name, d) in data["drag"])
#     waves = WaveParameters(; data["waves"]...)
#     materials = Dict(name => Material(d["type"], d["properties"]) for (name, d) in data["materials"])
#     motions = Dict(id => Motion(d["type"], d["properties"]) for (id, d) in data["motions"])
#     seabed = SeaBedParameters(; data["seabed"]...)

#     return ParameterHandler(points, segments, lines, drag, waves, materials, motions, seabed)
# end

# function to_dict(ph::ParameterHandler)::Dict
#     return Dict(
#         "points"    => Dict(id => Dict("tag" => p.tag, "coordinates" => p.coordinates,
#                                        "motion" => p.motion, "mesh_size" => p.mesh_size) for (id,p) in ph.points),
#         "segments"  => Dict(id => Dict("tag"=>s.tag, "start"=>s.start_point,
#                                        "stop"=>s.stop_point, "length"=>s.length,
#                                        "material"=>s.material, "density"=>s.density, "area"=>s.area) for (id,s) in ph.segments),
#         "lines"     => Dict(id => Dict("points"=>l.points, "segments"=>l.segments) for (id,l) in ph.lines),
#         "drag"      => Dict(name => Dict(field=>getfield(d, field) for field in fieldnames(DragProperties)) for (name,d) in ph.drag),
#         "waves"     => Dict(field=>getfield(ph.waves, field) for field in fieldnames(WaveParameters)),
#         "materials" => Dict(name => Dict("type"=>m.type, "properties"=>m.properties) for (name,m) in ph.materials),
#         "motions"   => Dict(id => Dict("type"=>m.type, "properties"=>m.properties) for (id,m) in ph.motions),
#         "seabed"    => Dict(field=>getfield(ph.seabed, field) for field in fieldnames(SeaBedParameters))
#     )
# end

# # --- Query/Modify ---
# function get_parameter(ph::ParameterHandler, category::Symbol, key)
#     return getfield(ph, category)[key]
# end

# function set_parameter!(ph::ParameterHandler, category::Symbol, key, value)
#     getfield(ph, category)[key] = value
# end

end # module
