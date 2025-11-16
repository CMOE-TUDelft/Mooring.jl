module ParameterHandlers

using Parameters
using YAML
using JSON3

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
  tag::String = "Point_$id"
  coords::Vector{Float64} = [0.0, 0.0]
  motion_tag::String = "default_motion"
  mesh_size::Float64 = 1.0
end

"""
PointParameters constructor

This constructor creates a new instance of the PointParameters struct with the given ID and optional parameters.
"""
function PointParameters(id::Int; coords::Vector{Float64}=[0.0,0.0], motion_tag::String="default", mesh_size::Float64=1.0)
  tag = "Point_$id"
  return PointParameters(id, tag, coords, motion_tag, mesh_size)
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
  tag::String = "Segment_$id"
  start_point::Int = 1
  stop_point::Int = 2
  length::Float64 = 10.0
  density::Float64 = 7850.0
  area::Float64 = 0.01
  material_tag::String = "default_material"
  drag_tag::String = "default_drag"
  seabed_tag::String = "default_seabed"
end

"""
SegmentParameters constructor

This constructor creates a new instance of the SegmentParameters struct with the given ID and optional parameters.
"""
function SegmentParameters(id::Int; start_point::Int=1, stop_point::Int=2, length::Float64=10.0,
  density::Float64=7850.0, area::Float64=0.01, material_tag::String="default_material",
  drag_tag::String="default_drag", seabed_tag::String="default_seabed")
  tag = "Segment_$id"
  SegmentParameters(id, tag, start_point, stop_point, length, density, area, material_tag, drag_tag, seabed_tag)
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
@with_kw struct LineParameters
  id::Int = 1
  tag::String = "Line_$id"
  points::Vector{Int} = [1,2]
  segments::Vector{Int} = [1]
end

"""
LineParameters constructor

This constructor creates a new instance of the LineParameters struct with the given ID and optional parameters.
"""
function LineParameters(id::Int; points::Vector{Int}=[1,2], segments::Vector{Int}=[1])
  tag = "Line_$id"
  LineParameters(id, tag, points, segments)
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
  f::String = "(t,x) -> VectorValue(0.0, 0.0)" # function f(t,x,y,z) of time (t) and position (x,y,z) in string format
  # WaveMotion parameters
  wave_tag::String = "default_waves"
end

"""
SeaBedParameters Struct
  
This struct contains the properties of the seabed.
The following parameters are included, with default values:
- tag::String = "default_seabed": Identifier for the seabed
- `kn::Real = 30e3`: Normal stiffness [N/m2]
- `linear_damping_factor::Real = 0.05`: Linear damping ratio [s]
- `quadratic_damping_factor::Real = 0.0`: Quadratic damping ratio [s^2/m]
- `od::Real = 0.1`: Outer diameter of the line [m]
- `A::Real = 0.008`: Area of the line [m^2]
- `tanh_ramp::Real = 1e2`: Tanh ramp function parameter 
- `penetration_depth_ramp::Real = 1e-3`: Penetration depth ramp function parameter [m]
- `still_weight::Real = 0.0`: Still weight [N]
- `cnstz::Real = 0.0`: Constant spring stiffness of the sea bed [N/m]
  
Relevant references:
- Quadratic law impact damping: https://doi.org/10.1080/0020739X.2021.1954253
- Critical damping of Moordyn: https://moordyn.readthedocs.io/en/latest/troubleshooting.html#model-stability-and-segment-damping
"""
@with_kw struct SeaBedParameters
  tag::String = "default_seabed"
  kn::Real = 30e3
  linear_damping_factor::Real = 0.05
  quadratic_damping_factor::Real = 0.0  
  od::Real = 0.1
  A::Real = 0.008 
  tanh_ramp::Real = 1.0e2
  penetration_depth_ramp::Real = 1.0e-3
  still_weight::Real = 0.0
  cnstz::Real = kn * od / A
end

"""
ParameterHandler
  
Central container for all parameters in a mooring system experiment.
It aggregates all other parameter structs into a single unified object.
  
# Fields
- `points::Dict{Int, PointParameters}` : Dictionary of point parameters keyed by point ID
- `segments::Dict{Int, SegmentParameters}` : Dictionary of segment parameters keyed by segment ID
- `lines::Dict{Int, LineParameters}` : Dictionary of line parameters keyed by line ID
- `drags::Dict{String, DragParameters}` : Dictionary of drag properties keyed by drag tag
- `waves::Dict{String, WaveParameters}` : Dictionary of wave parameter sets keyed by wave tag
- `materials::Dict{String, MaterialParameters}` : Dictionary of material models keyed by material tag
- `motions::Dict{String, MotionParameters}` : Dictionary of motions keyed by motion tag
- `seabeds::Dict{String, SeaBedParameters}` : Dictionary of seabed parameter sets keyed by tag
  
# Usage
```julia
ph = ParameterHandler()
  
# Add a point
ph.points[1] = PointParameters(id=1, coords=[0.0,0.0])
  
# Add seabed
ph.seabeds["default"] = SeaBedParameters()
  
# Query material
steel = ph.materials["steel"]
"""
mutable struct ParameterHandler
  points::Dict{Int, PointParameters}
  segments::Dict{Int, SegmentParameters}
  lines::Dict{Int, LineParameters}
  drags::Dict{String, DragParameters}
  waves::Dict{String, WaveParameters}
  materials::Dict{String, MaterialParameters}
  motions::Dict{String, MotionParameters}
  seabeds::Dict{String, Union{SeaBedParameters, Nothing}}
end

"""
ParameterHandler constructor

Create a ParameterHandler with default parameters.
"""
function ParameterHandler()
  return ParameterHandler(
  Dict{Int, PointParameters}(1=>PointParameters()),
  Dict{Int, SegmentParameters}(1=>SegmentParameters()),
  Dict{Int, LineParameters}(1=>LineParameters()),
  Dict{String, DragParameters}("default_drag"=>DragParameters()),
  Dict{String, WaveParameters}("default_waves"=>WaveParameters()),
  Dict{String, MaterialParameters}("default_material"=>MaterialParameters()),
  Dict{String, MotionParameters}("default_motion"=>MotionParameters()),
  Dict{String, Union{SeaBedParameters, Nothing}}("default_seabed"=>nothing)
  )
end

# -------------------------------
# YAML I/O
# -------------------------------

"""
    load_from_yaml(path::String) -> ParameterHandler
  
Load a YAML file defining experiment parameters into a `ParameterHandler`.
"""
function load_from_yaml(path::String)
  data = YAML.load_file(path)
  return _dict_to_handler(data)
end

"""
    save_to_yaml(ph::ParameterHandler, path::String)
  
Save the current parameter handler into a YAML file.
"""
function save_to_yaml(ph::ParameterHandler, path::String)
  open(path, "w") do io
    YAML.write(io, _handler_to_dict(ph))
  end
end

# -------------------------------
# JSON I/O
# -------------------------------

"""
    load_from_json(path::String) -> ParameterHandler
  
Load a JSON file defining experiment parameters into a `ParameterHandler`.
"""
function load_from_json(path::String)
    raw = read(path, String)
    data = JSON3.read(raw)
    dict_data = _json_to_dict(data)   # recursive conversion
    return _dict_to_handler(dict_data)
end

"""
    save_to_json(ph::ParameterHandler, path::String)
  
Save the current parameter handler into a JSON file.
"""
function save_to_json(ph::ParameterHandler, path::String)
  # json_data = _dict_to_json(_handler_to_dict(ph))
  json_data = _handler_to_dict(ph)
  json_str = JSON3.write(json_data; indent=4)
  open(path, "w") do io
    write(io, json_str)
  end
end

"""
    _json_to_dict(x)

Recursively convert `JSON3.Object` / `JSON3.Array` to plain Julia
`Dict{String,Any}` / `Vector{Any}` for uniform handling.
"""
function _json_to_dict(x)
    if x isa JSON3.Object
        return Dict(string(k) => _json_to_dict(v) for (k,v) in pairs(x))
    elseif x isa JSON3.Array
        return [_json_to_dict(v) for v in x]
    else
        return x
    end
end

"""
  _dict_to_json(x)

Recursively convert a Julia dictionary or array to a JSON3-compatible format.
"""
function _dict_to_json(x)
    if x isa Dict
        return JSON3.Object(Dict(string(k) => _dict_to_json(v) for (k,v) in pairs(x)))
    else
        return x
    end
end

# -------------------------------
# Internal converters
# -------------------------------

"""
    _dict_to_handler(data::Dict) -> ParameterHandler
  
Convert a parsed YAML/JSON dictionary into a `ParameterHandler`.
"""
function _dict_to_handler(data::Dict)
  ph = ParameterHandler()
    
  # Points
  if haskey(data, "points")
    for p in data["points"]
      pp = PointParameters(;_dict2kwargs(p)...)
      ph.points[pp.id] = pp
    end
  end
  
  # Segments
  if haskey(data, "segments")
    for s in data["segments"]
      sp = SegmentParameters(;_dict2kwargs(s)...)
      ph.segments[sp.id] = sp
    end
  end
  
  # Lines
  if haskey(data, "lines")
    for l in data["lines"]
      lp = LineParameters(; _dict2kwargs(l)...)
      ph.lines[lp.id] = lp
    end
  end
  
  # Drags
  if haskey(data, "drags")
    for d in data["drags"]
      dp = DragParameters(; _dict2kwargs(d)...)
      ph.drags[dp.tag] = dp
    end
  end
  
  # Waves
  if haskey(data, "waves")
    for w in data["waves"]
      wp = WaveParameters(; _dict2kwargs(w)...)
      ph.waves[wp.tag] = wp
    end
  end
  
  # Materials
  if haskey(data, "materials")
    for m in data["materials"]
      mp = MaterialParameters(; _dict2kwargs(m)...)
      ph.materials[mp.tag] = mp
    end
  end
  
  # Motions
  if haskey(data, "motions")
    for m in data["motions"]
      mo = MotionParameters(; _dict2kwargs(m)...)
      ph.motions[mo.tag] = mo
    end
  end
  
  # Seabed
  if haskey(data, "seabeds")
    for sb in data["seabeds"]
      sbo = SeaBedParameters(; _dict2kwargs(sb)...)
      ph.seabeds[sb["tag"]] = sbo
    end
  end
  
  return ph
end


"""
    _handler_to_dict(ph::ParameterHandler) -> Dict
  
Convert a `ParameterHandler` back into a dictionary for YAML/JSON export.
"""
function _handler_to_dict(ph::ParameterHandler)
  return Dict(
  "points"    => [Dict(String(field) => getfield(p, field) for field in fieldnames(PointParameters)) for p in values(ph.points)],
  "segments"  => [Dict(String(field) => getfield(s, field) for field in fieldnames(SegmentParameters)) for s in values(ph.segments)],
  "lines"     => [Dict(String(field) => getfield(l, field) for field in fieldnames(LineParameters)) for l in values(ph.lines)],
  "drags"     => [Dict(String(field) => getfield(d, field) for field in fieldnames(DragParameters)) for d in values(ph.drags)],
  "waves"     => [Dict(String(field) => getfield(w, field) for field in fieldnames(WaveParameters)) for w in values(ph.waves)],
  "materials" => [Dict(String(field) => getfield(m, field) for field in fieldnames(MaterialParameters)) for m in values(ph.materials)],
  "motions"   => [Dict(String(field) => getfield(mo, field) for field in fieldnames(MotionParameters)) for mo in values(ph.motions)],
  "seabeds"   => [Dict(String(field) => getfield(sb, field) for field in fieldnames(SeaBedParameters)) for sb in values(ph.seabeds) if sb !== nothing],
  )
end

function _dict2kwargs(d::AbstractDict)
    return [(Symbol(string(k)) => v) for (k,v) in pairs(d)]
end

end # module
