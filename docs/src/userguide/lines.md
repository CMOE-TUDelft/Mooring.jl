# Defining Lines

Mooring lines are **assemblies of mooring segments and points** that define the geometry and behavior of a complete line in the mooring system.  
Each line connects (multiple) anchors to (multiple) fairleads, or intermediate points, through one or more segments.

Each line is defined by:
- `id::Ing` - unique integer identifier
- `tag::String` – unique identifier of the line
- `points::Vector{Int}` – list of point IDs that belong to this line  
- `segments::Vector{Int}` – list of segment IDs that form the line  

---

## YAML example

```yaml
lines:
  - id: 1
    tag: line1
    points: [1, 2]
    segments: [1]
```
In this example:

A single line (`line1`) is created. It uses two points ("anchor" at the seabed and "fairlead" at the floater). These are connected by one segment (`seg1`, see example in [Mooring Segments](segments.md)).

---

## Julia example

```julia
using Mooring
import Mooring.ParameterHandler as PH
ph.points[1] = PH.PointParameters(id=1, coords=[0.0,0.0,-50.0], motion_tag="fixed")
ph.points[2] = PH.PointParameters(id=2, coords=[10.0,0.0,0.0], motion_tag="fairlead_motion")
ph.segments[1] = PH.SegmentParameters(id=1, start_point=1, stop_point=2, material_tag="steel", area=0.01,
                                      density=7850.0, length=10.0, drag_tag="no_drag", seabed_tag="default_seabed")
ph.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
```

---

Lines may consist of several segments joined together. Example:

```yaml
points:
  - id: 1
    tag: anchor
    coords: [0.0, 0.0, -100.0]
  - id: 2
    tag: midpoint
    coords: [25.0, 0.0, -50.0]
  - id: 3
    tag: fairlead
    coords: [50.0, 0.0, 0.0]

segments:
  - id: 1
    tag: lower
    start_point: 1
    stop_point: 2
    material: chain
    area: 0.007
    density: 7700.0
  - id: 2
    tag: upper
    start_point: 2
    stop_point: 3
    material: polyester
    area: 0.004
    density: 1400.0
lines:
  - id: 1
    tag: hybrid_line
    points: [1, 2, 3]
    segments: [1, 2]
```
This defines a hybrid line:

- Lower segment = chain (anchor → midpoint).
- Upper segment = polyester rope (midpoint → fairlead).

Lines are the system-level building blocks and are used to define [`MooringLine`](@ref) objects, which are the highest level entities in the simulation workflow.