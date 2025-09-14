# Defining Segments

Mooring segments represent the physical elements (ropes, chains, or synthetic cables) that connect two mooring points. Each segment connects a **start point** and a **stop point**, and has associated material and hydrodynamic properties.

Each segment is defined by the following parameters, see ([`SegmentParameters](@ref)):
- `id::Int` - unique integer identifier
- `tag::String` – unique name or identifier
- `start_point::Int` – ID of the starting point (refers to `points`)
- `stop_point::Int` – ID of the ending point (refers to `points`)
- `material_tag::String` – reference to a material tag (defined in `materials`)
- `area::Float64` – cross-sectional area [m²]
- `density::Float64` – mass density
- (optional) `drag_tag::String` – reference to a drag law (defined in `drags`)
- (optional) `seabed_tag::String` – reference to seabed interaction (defined in `seabeds`)

---

## YAML example

```yaml
segments:
  - id: 1
    tag: seg1
    start_point: 1
    stop_point: 2
    material: steel
    area: 0.01
    density: 7850.0
    drag: no_drag
    seabed: default_seabed
```

In this example:

`seg1` connects point 1 (e.g. an anchor) to point 2 (e.g. a fairlead).

It uses the "steel" material defined in `materials` (see materials section). A cross-sectional area of 0.01 m² and density 7850 kg/m³ are provided. Hydrodynamic drag and seabed interaction are also specified.

--- 

## Julia example

The same example can be directly defined in a julia script:

```julia
using Mooring
import Mooring.ParameterHandler as PH
ph.points[1] = PH.PointParameters(id=1, coords=[0.0,0.0,-50.0], motion_tag="fixed")
ph.points[2] = PH.PointParameters(id=2, coords=[10.0,0.0,0.0], motion_tag="fairlead_motion")
ph.segments[1] = PH.SegmentParameters(id=1, start_point=1, stop_point=2, material_tag="steel", area=0.01,
                                      density=7850.0, length=10.0, drag_tag="no_drag", seabed_tag="default_seabed")
```
Segments are turned into [`MooringSegment`](@ref) objects during line setup.