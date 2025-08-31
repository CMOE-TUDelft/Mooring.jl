module ParameterHandler

using Parameters

# using YAML
# using JSON


"""
 PointParameters

  This struct is used to define a point in the mooring system. It contains all relevant information
  about the point's position, motion, and mesh size. The following parameters are included, with default
  values provided for convenience:
    - `id::Int=1`: Unique identifier for the point.
    - `tag::String="Point_1"`: Human-readable name for the point.
    - `coords::Vector{Float64}=[0.0, 0.0]`: 2D coordinates of the point.
    - `motion_tag::String="default"`: Tag for the point's motion characteristics.
    - `mesh_size::Float64=1.0`: Mesh size for the point.
"""
@with_kw struct PointParameters
    id::Int = 1
    tag::String = "Point_1"
    coords::Vector{Float64} = [0.0, 0.0]
    motion_tag::String = "default_motion"
    mesh_size::Float64 = 1.0
    function PointParameters(id::Int; coords::Vector{Float64}=[0.0,0.0], motion_tag::String="default", mesh_size::Float64=1.0)
      tag = "Point_$id"
      new(id, tag, coords, motion_tag, mesh_size)
    end
end

"""
  SegmentParameters

This struct is used to define a segment in the mooring system. It contains all relevant information
about the segment properties. The following parameters and default values are used:
  - `id::Int=1`: Unique identifier for the segment.
  - `tag::String="Segment_1"`: Human-readable name for the segment.
  - `start_point::Int=1`: ID of the starting point.
  - `stop_point::Int=2`: ID of the stopping point.
  - `length::Float64=10.0`: Length of the segment.
  - `material::String="steel"`: Material of the segment.
  - `density::Float64=7850.0`: Density of the segment material.
  - `area::Float64=0.01`: Cross-sectional area of the segment.
  - `drag_tag::String="default"`: Tag for the segment's drag characteristics.
"""
@with_kw struct SegmentParameters
  id::Int = 1
  tag::String = "Segment_1"
  start_point::Int = 1
  stop_point::Int = 2
  length::Float64 = 10.0
  density::Float64 = 7850.0
  area::Float64 = 0.01
  material_tag::String = "default_material"
  drag_tag::String = "default_drag"
  function SegmentParameters(id::Int; start_point::Int=1, stop_point::Int=2, length::Float64=10.0,
    density::Float64=7850.0, area::Float64=0.01, material_tag::String="default_material", 
    drag_tag::String="default_drag")
    tag = "Segment_$id"
    new(id, tag, start_point, stop_point, length, density, area, material_tag, drag_tag)
  end
end

"""
LineParameters

  This struct is used to define a line in the mooring system. It contains the list of points and segments
  of the line. The following default values are used:
    - `id::Int=1`: Unique identifier for the line.
    - `tag::String="Line_1"`: Human-readable name for the line.
    - `points::Vector{Int}=[1,2]`: List of point IDs.
    - `segments::Vector{Int}=[1]`: List of segment IDs.
"""
struct LineParameters
  id::Int = 1
  tag::String = "Line_1"
  points::Vector{Int} = [1,2]
  segments::Vector{Int} = [1]
  function LineParameters(id::Int; points::Vector{Int}=[1,2], segments::Vector{Int}=[1])
    tag = "Line_$id"
    new(id, tag, points, segments)
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
@with_kw struct DragParameters
  tag::String = "default_drag"
  dragType::String = "NoDrag"
  ρw::Float64 = 1025.0
  nd::Float64 = 0.0
  od::Float64 = 0.0
  id::Float64 = 0.0
  AStr::Float64 = 0.0
  Cd_n::Float64 = 0.0
  Cd_t::Float64 = 0.0
  dd_n::Float64 = 0.0
  dd_t::Float64 = 0.0
  Ca_n::Float64 = 0.0
  Ca_t::Float64 = 0.0
  Cfd_n::Float64 = 0.0
  Cfd_t::Float64 = 0.0
  function DragParameters(; tag="default", dragType="NoDrag", ρw=0.0, nd=0.0, od=0.0, id=0.0, AStr=0.0,
    Cd_n=0.0, Cd_t=0.0, dd_n=0.0, dd_t=0.0,
    Ca_n=0.0, Ca_t=0.0, Cfd_n=0.0, Cfd_t=0.0)
    new(tag, dragType, ρw, nd, od, id, AStr, Cd_n, Cd_t, dd_n, dd_t, Ca_n, Ca_t, Cfd_n, Cfd_t)
  end
end

"""
WaveParameters

This struct contains the parameters for the wave conditions.
The following parameters are included, with default values:
- `Hs::Real = 0.0`: Significant wave height
- `Tp::Real = 0.0`: Peak wave period
- `h0::Real = 100.0`: Water depth
- `nω::Int = 64`: Number of frequency components
- `seed::Int = 0`: Seed for random phase
- `ωc::Real = -1.0`: Cut-off frequency
- `enableWaveSpec::Bool = false`: Enable wave spectrum
"""
@with_kw struct WaveParameters
  tag::String = "default_waves"
  Hs::Real = 0.0
  Tp::Real = 0.0
  h0::Real = 100.0
  nω::Int = 64
  seed::Int = 0
  ωc::Real = -1.0
  enableWaveSpec::Bool = false

end

"""
MaterialParameters

This struct contains the material properties for a given mooring segment. There are two implemented
material models: "LinearElastic" and "Scharpery". The default choice is "LinearElastic". The following
parameters and default values are used:
- tag::String = "default_material": Identifier for the material
- type::String = "LinearElastic": Material model type
For LinearElastic:
- E::Real = 1.0: Young's modulus
- ν::Real = 0.3: Poisson's ratio
For Scharpery:
- D0::Real = 1.0: Elastic compliance
- N::Int = 1: Number of relaxation times
- Dn::Vector{Real} = [1.0]: Compliance for each relaxation time
- λn::Vector{Real} = [1.0]: Inverse of relaxation times
- g0::String = "1.0": function of σ in string format
- g1::String = "1.0": function of σ in string format
- g2::String = "1.0": function of σ in string format
"""
@with_kw struct MaterialParameters
  tag::String = "default_material"
  type::String = "LinearElastic"
  # Linear elastic properties
  E::Real = 1.0
  μ::Real = 0.5*E # second Lame constant
  # Scharpery properties
  D0::Real = 1.0
  N::Int = 1
  Dn::Vector{Real} = [1.0]
  λn::Vector{Real} = [1.0]
  g0::String = "1.0"
  g1::String = "1.0"
  g2::String = "1.0"
end

"""
MotionParameters

This struct contains the parameters required to describe a point motion. The following motion 
types are supported:
- `CustomMotion`: User-defined motion with custom function.
- `WaveMotion`: Motion driven by wave parameters.
- `Nothing`: No motion applied, the point is free to move.
The following parameters and default values are used:
- tag::String = "default_motion": Identifier for the motion
- type::String = "CustomMotion": Type of motion
for `CustomMotion`:
- f::String = "0.0": function f(t,x,y,z) of time (t) and position (x,y,z) in string format
for `WaveMotion`:
- wave_tag::String = "default_waves": Identifier for the wave motion
"""
@with_kw struct MotionParameters
  tag::String = "default_motion"
  type::String = "CustomMotion"
  # CustomMotion parameters
  f::String = "0.0" # function f(t,x,y,z) of time (t) and position (x,y,z) in string format
  # WaveMotion parameters
  wave_tag::String = "default_waves"
end

#### HERE!


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
