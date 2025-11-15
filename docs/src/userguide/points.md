# Defining Points

Points are fundamental entities in **Mooring.jl** that define the endpoints of segments. They represent various connection points in a mooring system:
- **Anchors**: Fixed points on the seabed
- **Fairleads**: Connection points on floating structures
- **Intermediate connections**: Joints between different mooring components

#### Table of Contents

- [Two Types of Point Structures](#two-types-of-point-structures)
- [Point Properties (PointParameters)](#point-properties-pointparameters)
- [Defining Points directly in Julia script](#defining-points-directly-in-julia-script)
- [Defining Points in YAML input file](#defining-points-in-yaml-input-file)
- [Working with Points](#working-with-points)
- [Workflow: From PointParameters to MooringPoint](#workflow-from-pointparameters-to-mooringpoint)
- [Common Issues](#common-issues)
- [Next Steps](#next-steps)
- [See Also](#see-also)

## Two Types of Point Structures

**Mooring.jl** uses two distinct structures for handling points at different levels:

### 1. `PointParameters` (High-Level: Parameter Handling)

`PointParameters` is used for **defining and configuring** points in YAML files or through the parameter handler. This struct stores all the user-specified properties needed to create a point in the mooring system. You interact with this struct when setting up your simulation configuration.

```julia
using Mooring.ParameterHandlers as PH

# Create point parameters for configuration
point_params = PH.PointParameters(
    id=1,
    coords=[0.0, -100.0],
    motion_tag="fixed",
    mesh_size=10.0
)
```

### 2. `MooringPoint` (Low-Level: FEM Operations)

`MooringPoint` is used internally for **finite element analysis** operations. This struct contains the boundary triangulation and motion functions needed for FEM computations. The package automatically creates `MooringPoint` instances from your `PointParameters` when building the discrete model.

```julia
using Mooring.MooringPoints

# MooringPoint is created internally from the discrete model
# You typically don't create these directly
mooring_point = MooringPoint(model, "Point_1", motion_type)
```

**Key Distinction:**
- **`PointParameters`**: What you define (coordinates, properties) → *User configuration level*
- **`MooringPoint`**: What the solver uses (triangulation, motion functions) → *FEM computation level*

## Point Properties (PointParameters)

When defining points through the parameter handler, each point uses the following parameters:

| Parameter | Type | Description |
|-----------|------|-------------|
| `id` | `Int` | Unique integer identifier for the point |
| `coords` | `Vector{Float64}` | Spatial coordinates in meters: `[x, z]` for 2D or `[x, y, z]` for 3D |
| `motion_tag` | `String` | Reference to a defined motion type |
| `mesh_size` | `Float64` | Controls local mesh refinement (optional) |

### Coordinate System

Mooring.jl supports both **2D** and **3D** coordinate systems:

#### 2D Coordinates
For planar mooring problems, use `[x, z]` coordinates:
- **x**: Horizontal coordinate [m]
- **z**: Vertical coordinate [m] (positive upward)
  - z = 0: Water surface (typical)
  - z < 0: Below water surface
  - z = -depth: Seabed level

```julia
# 2D example: vertical plane analysis
ph.points[1] = PH.PointParameters(
    id=1, 
    coords=[0.0, -100.0],  # x=0m, z=-100m (seabed)
    motion_tag="fixed", 
    mesh_size=10.0
)
```

#### 3D Coordinates
For full spatial mooring systems, use `[x, y, z]` coordinates:
- **x, y**: Horizontal plane coordinates [m]
- **z**: Vertical coordinate [m] (positive upward)
  - z = 0: Water surface (typical)
  - z < 0: Below water surface
  - z = -depth: Seabed level

```julia
# 3D example: multi-directional mooring
ph.points[1] = PH.PointParameters(
    id=1, 
    coords=[50.0, 30.0, -100.0],  # x=50m, y=30m, z=-100m
    motion_tag="fixed", 
    mesh_size=10.0
)
```

### Motion Types

Points can have different motion constraints:

1. **Fixed**: Point is stationary in space
   ```julia
   motion_tag = "fixed"
   ```

2. **Prescribed Motion**: Point follows a user-defined trajectory
   ```julia
   motion_tag = "vessel_motion"  # Custom motion function
   ```

3. **Free**: Point can move (for coupled analysis)
   ```julia
   motion_tag = "free"  # Future feature
   ```

### Mesh Size

The `mesh_size` parameter controls the finite element mesh density near the point:
- **Smaller values** (e.g., 1-5 m): Fine mesh, higher accuracy, slower
- **Larger values** (e.g., 10-20 m): Coarse mesh, lower accuracy, faster
- **Recommendation**: Use finer mesh where you expect large gradients (e.g., near fairleads)

---

## Defining Points directly in Julia script

### Basic 2D Example

For planar analysis (e.g., single mooring line in vertical plane):

```julia
using Mooring
import Mooring.ParameterHandlers as PH

# Create parameter handler
ph = PH.ParameterHandler()

# Define motion types
ph.motions["fixed"] = PH.MotionParameters(tag="fixed")

# Define an anchor point (2D)
ph.points[1] = PH.PointParameters(
    id=1, 
    coords=[0.0, -100.0],  # x=0m, z=-100m below surface
    motion_tag="fixed", 
    mesh_size=10.0
)

# Define a fairlead point (2D)
ph.points[2] = PH.PointParameters(
    id=2, 
    coords=[50.0, -20.0],  # x=50m horizontal, z=-20m below surface
    motion_tag="fixed", 
    mesh_size=5.0  # Finer mesh near fairlead
)
```

### Basic 3D Example

For full spatial analysis:

```julia
using Mooring
import Mooring.ParameterHandlers as PH

# Create parameter handler
ph = PH.ParameterHandler()

# Define motion types
ph.motions["fixed"] = PH.MotionParameters(tag="fixed")

# Define an anchor point (3D)
ph.points[1] = PH.PointParameters(
    id=1, 
    coords=[0.0, 0.0, -100.0],  # x=0m, y=0m, z=-100m below surface
    motion_tag="fixed", 
    mesh_size=10.0
)

# Define a fairlead point (3D)
ph.points[2] = PH.PointParameters(
    id=2, 
    coords=[50.0, 30.0, -20.0],  # x=50m, y=30m, z=-20m below surface
    motion_tag="fixed", 
    mesh_size=5.0  # Finer mesh near fairlead
)
```

### Three-Line Mooring System

Example of a symmetric three-line mooring configuration:

```julia
using Mooring
import Mooring.ParameterHandlers as PH

ph = PH.ParameterHandler()
ph.motions["fixed"] = PH.MotionParameters(tag="fixed")

# System parameters
depth = 600.0                    # Water depth [m]
angles = [π/3, π, 5π/3]         # 60°, 180°, 300°
rAnchor = 1600.0                # Anchor radius [m]
zFair = -21.0                   # Fairlead depth [m]
rFair = 20.0                    # Fairlead radius [m]

# Central point (optional - for visualization)
ph.points[1] = PH.PointParameters(
    id=1,
    coords=[0.0, 0.0, zFair],
    motion_tag="fixed",
    mesh_size=10.0
)

point_id = 2
# Create anchor and fairlead points for each line
for (i, angle) in enumerate(angles)
    # Anchor points
    anchor_x = rAnchor * cos(angle)
    anchor_y = rAnchor * sin(angle)
    ph.points[point_id] = PH.PointParameters(
        id=point_id,
        coords=[anchor_x, anchor_y, -depth],
        motion_tag="fixed",
        mesh_size=20.0
    )
    point_id += 1
    
    # Fairlead points
    fair_x = rFair * cos(angle)
    fair_y = rFair * sin(angle)
    ph.points[point_id] = PH.PointParameters(
        id=point_id,
        coords=[fair_x, fair_y, zFair],
        motion_tag="fixed",
        mesh_size=5.0
    )
    point_id += 1
end
```

---

## Defining Points in YAML input file

For larger systems, YAML configuration files are more convenient:

### Simple 2D System

For planar mooring analysis:

```yaml
# mooring_config_2d.yaml
points:
  - id: 1
    coords: [0.0, -100.0]      # 2D: x, z coordinates
    motion_tag: "fixed"
    mesh_size: 10.0
    
  - id: 2
    coords: [50.0, -20.0]      # 2D: x, z coordinates
    motion_tag: "fixed"
    mesh_size: 5.0
```

### Simple 3D System

For full spatial analysis:

```yaml
# mooring_config_3d.yaml
points:
  - id: 1
    coords: [0.0, 0.0, -100.0]   # 3D: x, y, z coordinates
    motion_tag: "fixed"
    mesh_size: 10.0
    
  - id: 2
    coords: [50.0, 30.0, -20.0]  # 3D: x, y, z coordinates
    motion_tag: "fixed"
    mesh_size: 5.0
```

Load in Julia:
```julia
using Mooring
import Mooring.ParameterHandlers as PH

ph = PH.read_parameters("mooring_config.yaml")
```

### Multi-Line Configuration

```yaml
points:
  # Anchors
  - id: 1
    coords: [800.0, 1385.6, -600.0]  # 60° angle
    motion_tag: "fixed"
    mesh_size: 20.0
    
  - id: 2
    coords: [-1600.0, 0.0, -600.0]   # 180° angle
    motion_tag: "fixed"
    mesh_size: 20.0
    
  - id: 3
    coords: [800.0, -1385.6, -600.0] # 300° angle
    motion_tag: "fixed"
    mesh_size: 20.0
    
  # Fairleads
  - id: 4
    coords: [10.0, 17.32, -21.0]     # 60° angle
    motion_tag: "vessel_motion"
    mesh_size: 5.0
    
  - id: 5
    coords: [-20.0, 0.0, -21.0]      # 180° angle
    motion_tag: "vessel_motion"
    mesh_size: 5.0
    
  - id: 6
    coords: [10.0, -17.32, -21.0]    # 300° angle
    motion_tag: "vessel_motion"
    mesh_size: 5.0

motions:
  - tag: "fixed"
    type: "NoMotion"
    
  - tag: "vessel_motion"
    type: "CustomMotion"
    wave_tag: "(t,(x,y))->VectorValue(0.0, 0.0)"
```

---

## Working with Points

### Accessing Point Data

```julia
# Get a specific point
point_1 = ph.points[1]

# Access properties
coords = point_1.coords
motion = point_1.motion_tag
mesh = point_1.mesh_size

println("Point 1 at: $coords")
```

### Iterating Over Points

```julia
# Loop through all points
for (id, point) in ph.points
    println("Point $id: $(point.coords)")
end

# Count points
num_points = length(ph.points)
println("Total points: $num_points")
```

### Modifying Points

```julia
# Update point coordinates (e.g., adjust fairlead position)
ph.points[2] = PH.PointParameters(
    id=2,
    coords=[55.0, 0.0, -22.0],  # New position
    motion_tag="fixed",
    mesh_size=5.0
)

# Or modify in place
old_point = ph.points[2]
ph.points[2] = PH.PointParameters(
    id=old_point.id,
    coords=old_point.coords .+ [5.0, 0.0, -2.0],  # Shift position
    motion_tag=old_point.motion_tag,
    mesh_size=old_point.mesh_size
)
```

---

## Workflow: From PointParameters to MooringPoint

Understanding how `PointParameters` transforms into `MooringPoint` helps clarify the package architecture:

### 1. **User Configuration Stage** (`PointParameters`)

You define point properties in YAML or Julia:

```julia
# Define configuration
ph.points[1] = PH.PointParameters(
    id=1,
    coords=[0.0, -100.0],
    motion_tag="fixed",
    mesh_size=10.0
)
```

At this stage, you're working with **physical properties** (coordinates, motion type, mesh size).

### 2. **Geometry Generation Stage**

The package uses your `PointParameters` to create a discrete geometric model:

```julia
# Package internally creates geometry from your parameters
model = create_mooring_discrete_model(ph)
```

The coordinates and mesh sizes from `PointParameters` are used to generate the mesh geometry.

### 3. **FEM Setup Stage** (`MooringPoint`)

The package creates `MooringPoint` instances for FEM operations:

```julia
# Package internally creates MooringPoint from the discrete model
mooring_point = MooringPoint(
    model,              # Discrete model with mesh
    "Point_1",          # Tag from PointParameters
    motion_type         # Motion function from motion_tag
)
```

At this stage, you're working with **FEM structures** (triangulations, motion functions).

### 4. **What Each Structure Contains**

| Aspect | PointParameters | MooringPoint |
|--------|----------------|--------------|
| **Purpose** | User configuration | FEM computation |
| **Contains** | Coordinates, mesh size, motion tag | Boundary triangulation, motion functions |
| **When used** | Parameter setup, YAML I/O | Assembly, solving, post-processing |
| **Typical user** | You (simulation setup) | Package internals (solver) |
| **Mutability** | Can be modified | Typically fixed once created |

### 5. **When You Interact With Each**

**Work with `PointParameters` when:**
- Reading/writing YAML configuration files
- Setting up simulation parameters
- Modifying point locations before solving
- Parametric studies (changing coordinates, mesh sizes)

**Work with `MooringPoint` when:**
- Extending the package with custom motion types
- Implementing custom boundary conditions
- Developing advanced FEM features
- Debugging low-level FEM operations

**Key Takeaway:** You primarily work with `PointParameters` for configuration. The package automatically handles the conversion to `MooringPoint` for FEM computations.

---

## Common Issues

### Problem: Points too close together

**Symptoms**: Mesh generation fails or very small elements
```julia
# Bad: Points only 0.1m apart
ph.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0, 0.0], ...)
ph.points[2] = PH.PointParameters(id=2, coords=[0.1, 0.0, 0.0], ...)
```

**Solution**: Ensure minimum separation (> mesh_size)
```julia
# Good: Points 10m apart with 5m mesh
ph.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0, 0.0], mesh_size=5.0, ...)
ph.points[2] = PH.PointParameters(id=2, coords=[10.0, 0.0, 0.0], mesh_size=5.0, ...)
```

---

## Next Steps

After defining points:
1. Define [Segments](segments.md) connecting your points
2. Configure [Material Properties](../API/Physics/Materials.md)
3. Assemble [Complete Lines](lines.md)

## See Also

- [API: MooringPoint](../API/Entities/MooringPoint.md)
- [API: PointMotion](../API/Physics/PointMotion.md)
- [API: ParameterHandler](../API/IO/ParameterHandler.md)