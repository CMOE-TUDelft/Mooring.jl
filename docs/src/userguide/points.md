# Defining Points

Points define the ends of segments. They represent fixed or moving anchors, fairleads or connection points in a mooring system.

They can be fixed (e.g. anchors on the seabed) or moving (e.g. connection points on a floating platform).  

Each point is defined by:
- `id::Int` - a unique integer identifier
- `tag::String` – a unique name or identifier
- `coords::Vector{Float64}` – spatial coordinates `[x,y,z]` in meters
- `motion::String` – type name of the motion constraint (`"fixed"`, `"prescribed"`, or `"free"`)
- `mesh_size::Float` - (optional) additional parameters to prescribe mesh size

See [`PointParameters`](@ref) for a description of all possible parameters.

---

## YAML Example
We can define a point through the YAML parameter handler. Here is an example of two points named "anchor" and "fairlead", for which we specify the motion tags "fixed_motion" and "fairlead_motion". Note that these motions are defined through the [`PointMotion`](@ref) type, see also [*Point Motion*](userguide/motion.md).
```yaml
points:
  - id: 1
    tag: anchor
    coords: [0.0, 0.0, -50.0]
    motion_tag: "fixed"
  - id: 2
    tag: fairlead
    coords: [10.0, 0.0, 0.0]
    motion_tag: "fairlead_motion"
```

---

## Julia Example
The same points can be directly defined in a julia script.
```julia
using Mooring
import Mooring.ParameterHandler as PH
ph.points[1] = PH.PointParameters(id=1, coords=[0.0,0.0,-50.0], motion_tag="fixed")
ph.points[2] = PH.PointParameters(id=2, coords=[10.0,0.0,0.0], motion_tag="fairlead_motion")
```


Points are turned into [`MooringPoint`](@ref) objects during line setup.