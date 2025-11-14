# Defining Points

Points are fundamental entities in **Mooring.jl** that define the endpoints of segments. They represent various connection points in a mooring system:
- **Anchors**: Fixed points on the seabed
- **Fairleads**: Connection points on floating structures
- **Intermediate connections**: Joints between different mooring components

## Point Properties

Each point is characterized by the following parameters:

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

## Defining Points in Julia

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

## Defining Points in YAML

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