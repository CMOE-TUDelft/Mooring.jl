# Defining Segments

Mooring segments are the fundamental physical elements in **Mooring.jl** that represent cables connecting two points. Segments can model:
- **Steel chains**: Heavy, rigid mooring chains
- **Wire ropes**: Steel wire cables with high strength
- **Synthetic ropes**: Lightweight polyester or nylon cables
- **Hybrid systems**: Combinations of different materials

Each segment connects a **start point** to a **stop point** and carries material, geometric, and hydrodynamic properties.

#### Table of Contents

- [Two Types of Segment Structures](#two-types-of-segment-structures)
- [Segment Properties (SegmentParameters)](#segment-properties-segmentparameters)
- [Defining Segments in Julia](#defining-segments-in-julia)
- [Defining Segments in YAML](#defining-segments-in-yaml)
- [Working with Segments](#working-with-segments)
- [Workflow: From SegmentParameters to MooringSegment](#workflow-from-segmentparameters-to-mooringsegment)
- [Physical Map Functions](#physical-map-functions)
- [Best Practices](#best-practices)
- [Common Issues](#common-issues)
- [Next Steps](#next-steps)
- [See Also](#see-also)

## Two Types of Segment Structures

**Mooring.jl** uses two distinct structures for handling segments at different levels:

### 1. `SegmentParameters` (High-Level: Parameter Handling)

`SegmentParameters` is used for **defining and configuring** segments in YAML files or through the parameter handler. This struct stores all the user-specified properties needed to create a segment (cable) in the mooring system. You interact with this struct when setting up your simulation configuration.

```julia
using Mooring.ParameterHandlers as PH

# Create segment parameters for configuration
segment_params = PH.SegmentParameters(
    id=1,
    start_point=1,
    stop_point=2,
    length=100.0,
    area=0.00785,
    density=7850.0,
    material_tag="steel_chain"
)
```

### 2. `MooringSegment` (Low-Level: FEM Operations)

`MooringSegment` is used internally for **finite element analysis** operations. This struct contains the triangulation and material properties. From these parameters it can create **its own dedicated finite element space**. The package automatically creates `MooringSegment` instances from your `SegmentParameters` when building the discrete model.

```julia
using Mooring.MooringSegments

# MooringSegment is created internally from the discrete model
# You typically don't create these directly
mooring_segment = MooringSegment(
    model,              # Discrete model with mesh
    "Segment_1",        # Tag from SegmentParameters
    pointA,             # MooringPoint instance
    pointB,             # MooringPoint instance
    map_function,       # Geometric map (linear/quadratic)
    material,           # Material properties
    density,            # Submerged density
    area                # Cross-sectional area
)

# Key feature: Each MooringSegment gets its own FE space
U, V = get_transient_FESpaces(mooring_segment, order=1, dim=2)
```

**Key Distinction:**
- **`SegmentParameters`**: What you define (geometry, material, connections) → *User configuration level*
- **`MooringSegment`**: What the solver uses (triangulation, FE space, residuals) → *FEM computation level*

**Critical Feature:** Each `MooringSegment` has its **own independent finite element space**. This design allows:
- Different mesh refinements per segment
- Different interpolation orders per segment
- Flexible boundary conditions at segment endpoints
- Efficient assembly of multi-segment systems

## Segment Properties (SegmentParameters)

When defining segments through the parameter handler, each segment uses the following parameters:

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `id` | `Int` | Yes | Unique integer identifier for the segment |
| `start_point` | `Int` | Yes | ID of the starting point |
| `stop_point` | `Int` | Yes | ID of the ending point |
| `material_tag` | `String` | Yes | Reference to a defined material |
| `area` | `Float64` | Yes | Cross-sectional area [m²] |
| `density` | `Float64` | Yes | Material mass density [kg/m³] |
| `length` | `Float64` | Yes | Unstretched length of the segment [m] |
| `drag_tag` | `String` | No | Reference to drag model (default: "no_drag") |
| `seabed_tag` | `String` | No | Reference to seabed interaction (default: "default_seabed") |

### Unstretched Length

The `length` parameter specifies the **unstretched** (reference) length of the cable:
- This is the natural length without any applied tension
- Under tension, the actual length will increase based on material stiffness
- For catenary configurations, `length` should be greater than the straight-line distance between points

### Cross-Sectional Area

The area affects:
- **Axial stiffness**: Larger area → stiffer cable
- **Weight**: Weight per unit length = `density × area × g`
- **Hydrodynamic forces**: Drag is proportional to area (if drag enabled)

**Typical values**:
| Cable Type | Diameter | Area [m²] | Notes |
|------------|----------|-----------|-------|
| Chain 76mm | 76mm | 0.00454 | ≈ π×(0.038)² |
| Chain 120mm | 120mm | 0.0113 | ≈ π×(0.06)² |
| Wire rope 100mm | 100mm | 0.00785 | Circular cross-section |
| Synthetic rope 200mm | 200mm | 0.0314 | May be non-circular |

### Density

Density of the submerged material. The material density determines the submerged weight:
```
Submerged weight per meter = ρ_submerged × area × g = (ρ_cable - ρ_water) × area × g
```

**Typical densities**:
| Material | Density [kg/m³] | Notes |
|----------|-----------------|-------|
| Steel chain | 7850 | Standard structural steel |
| Wire rope | 7800-7900 | Slightly lower due to voids |
| Polyester rope | 1380 | Positive buoyancy in water |
| Nylon rope | 1140 | Positive buoyancy in water |
| Seawater | 1025 | For reference |

### Material Tag

References a material, which is defined depending on the material model. For an elastic material, with elastic properties:
- **Young's modulus (E)**: Axial stiffness
- **Shear modulus (μ)**: Lateral stiffness
- Supports linear and nonlinear material models

See [Materials](../API/Physics/Materials.md) for more details.

### Drag Tag

Optional reference to hydrodynamic drag model:
- `"no_drag"`: No hydrodynamic forces (default, for quasi-static in still water)
- Custom drag models: For dynamic analysis with currents/waves

See [Drag](../API/Physics/Drag.md) for more details.

### Seabed Tag

Optional reference to seabed contact model:
- `"default_seabed"`: Standard seabed contact model
- Custom models: For specific soil properties

See [SeaBed](../API/Physics/SeaBed.md) for more details.

---

## Defining Segments in Julia

### Basic 2D Example

Simple vertical hanging cable:

```julia
using Mooring
import Mooring.ParameterHandlers as PH

# Create parameter handler
ph = PH.ParameterHandler()

# Define material (steel chain)
ph.materials["chain"] = PH.MaterialParameters(
    tag="chain",
    E=2.1e11,    # Young's modulus [Pa]
    μ=8.1e10     # Shear modulus [Pa]
)

# Define motion type
ph.motions["fixed"] = PH.MotionParameters(tag="fixed")

# Define points (2D)
ph.points[1] = PH.PointParameters(
    id=1,
    coords=[0.0, -100.0],    # Anchor at seabed
    motion_tag="fixed",
    mesh_size=10.0
)
ph.points[2] = PH.PointParameters(
    id=2,
    coords=[0.0, -20.0],     # Fairlead
    motion_tag="fixed",
    mesh_size=5.0
)

# Define segment
ph.segments[1] = PH.SegmentParameters(
    id=1,
    start_point=1,
    stop_point=2,
    material_tag="chain",
    area=0.0113,              # 120mm chain ≈ 0.0113 m²
    density=7850.0,           # Steel density [kg/m³]
    length=90.0,              # Unstretched length [m]
    drag_tag="no_drag",
    seabed_tag="default_seabed"
)
```

### Catenary Cable (3D)

Cable with excess length for natural catenary shape:

```julia
using Mooring
import Mooring.ParameterHandlers as PH

ph = PH.ParameterHandler()

# Define material
ph.materials["wire_rope"] = PH.MaterialParameters(
    tag="wire_rope",
    E=1.8e11,    # Wire rope slightly lower than solid steel
    μ=7.0e10
)

ph.motions["fixed"] = PH.MotionParameters(tag="fixed")

# Define anchor point
ph.points[1] = PH.PointParameters(
    id=1,
    coords=[0.0, 0.0, -150.0],    # 150m depth
    motion_tag="fixed",
    mesh_size=15.0
)

# Define fairlead point
ph.points[2] = PH.PointParameters(
    id=2,
    coords=[100.0, 50.0, -25.0],  # 100m in x, 50m in y, 25m depth
    motion_tag="fixed",
    mesh_size=8.0
)

# Calculate distances
horizontal_dist = sqrt(100.0^2 + 50.0^2)  # ≈ 111.8 m
vertical_dist = abs(-25.0 - (-150.0))     # = 125 m
straight_dist = sqrt(horizontal_dist^2 + vertical_dist^2)  # ≈ 167.4 m

# Add 15% excess for catenary
segment_length = straight_dist * 1.15  # ≈ 192.5 m

# Define segment
ph.segments[1] = PH.SegmentParameters(
    id=1,
    start_point=1,
    stop_point=2,
    material_tag="wire_rope",
    area=0.00785,             # 100mm wire rope
    density=7850.0,
    length=segment_length,
    drag_tag="no_drag",
    seabed_tag="default_seabed"
)
```

### Multi-Segment Line

Chain-polyester-chain configuration:

```julia
using Mooring
import Mooring.ParameterHandlers as PH

ph = PH.ParameterHandler()

# Define materials
ph.materials["chain"] = PH.MaterialParameters(tag="chain", E=2.1e11, μ=8.1e10)
ph.materials["polyester"] = PH.MaterialParameters(tag="polyester", E=5.0e9, μ=2.0e9)

ph.motions["fixed"] = PH.MotionParameters(tag="fixed")

# Define points (3 points for 2 segments)
ph.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0, -200.0], motion_tag="fixed", mesh_size=20.0)  # Anchor
ph.points[2] = PH.PointParameters(id=2, coords=[100.0, 0.0, -150.0], motion_tag="fixed", mesh_size=10.0) # Intermediate
ph.points[3] = PH.PointParameters(id=3, coords=[150.0, 0.0, -30.0], motion_tag="fixed", mesh_size=8.0)   # Fairlead

# Bottom chain segment
ph.segments[1] = PH.SegmentParameters(
    id=1,
    start_point=1,
    stop_point=2,
    material_tag="chain",
    area=0.00785,             # 100mm chain
    density=7850.0,
    length=120.0,             # Ground chain
    drag_tag="no_drag",
    seabed_tag="default_seabed"
)

# Polyester rope segment
ph.segments[2] = PH.SegmentParameters(
    id=2,
    start_point=2,
    stop_point=3,
    material_tag="polyester",
    area=0.0314,              # 200mm polyester rope
    density=1380.0,           # Lighter than water!
    length=200.0,             # Main suspension
    drag_tag="no_drag",
    seabed_tag="default_seabed"
)
```

---

## Defining Segments in YAML

### Simple Single Segment

```yaml
# mooring_segment.yaml
materials:
  - tag: "chain"
    E: 2.1e11
    μ: 8.1e10

motions:
  - tag: "fixed"
    type: "NoMotion"

points:
  - id: 1
    coords: [0.0, -100.0]
    motion_tag: "fixed"
    mesh_size: 10.0
  - id: 2
    coords: [50.0, -20.0]
    motion_tag: "fixed"
    mesh_size: 5.0

segments:
  - id: 1
    start_point: 1
    stop_point: 2
    material_tag: "chain"
    area: 0.0113
    density: 7850.0
    length: 90.0
    drag_tag: "no_drag"
    seabed_tag: "default_seabed"
```

### Multi-Segment Configuration

```yaml
segments:
  # Bottom chain
  - id: 1
    start_point: 1
    stop_point: 2
    material_tag: "chain"
    area: 0.00785
    density: 7850.0
    length: 120.0
    drag_tag: "no_drag"
    seabed_tag: "default_seabed"
    
  # Polyester rope
  - id: 2
    start_point: 2
    stop_point: 3
    material_tag: "polyester"
    area: 0.0314
    density: 1380.0
    length: 200.0
    drag_tag: "no_drag"
    seabed_tag: "default_seabed"
    
  # Top chain
  - id: 3
    start_point: 3
    stop_point: 4
    material_tag: "chain"
    area: 0.00454
    density: 7850.0
    length: 50.0
    drag_tag: "no_drag"
    seabed_tag: "default_seabed"
```

---

## Working with Segments

### Accessing Segment Data

```julia
# Get a specific segment
segment_1 = ph.segments[1]

# Access properties
start_pt = segment_1.start_point
stop_pt = segment_1.stop_point
material = segment_1.material_tag
area = segment_1.area
length = segment_1.length

println("Segment 1: Point $start_pt → Point $stop_pt")
println("  Material: $material, Area: $area m², Length: $length m")
```

### Iterating Over Segments

```julia
# Loop through all segments
for (id, segment) in ph.segments
    weight_per_meter = segment.density * segment.area * 9.81
    println("Segment $id: $(weight_per_meter) N/m")
end

# Count segments
num_segments = length(ph.segments)
println("Total segments: $num_segments")
```

### Calculating Segment Properties

```julia
# Axial stiffness
function axial_stiffness(segment::PH.SegmentParameters, material::PH.MaterialParameters)
    return material.E * segment.area
end

# Example usage
seg = ph.segments[1]
mat = ph.materials[seg.material_tag]
EA = axial_stiffness(seg, mat)
println("Axial stiffness EA = $EA N")
```

---

## Workflow: From SegmentParameters to MooringSegment

Understanding how `SegmentParameters` transforms into `MooringSegment` helps clarify the package architecture and the role of independent FE spaces:

### 1. **User Configuration Stage** (`SegmentParameters`)

You define segment properties in YAML or Julia:

```julia
# Define configuration
ph.segments[1] = PH.SegmentParameters(
    id=1,
    start_point=1,
    stop_point=2,
    length=100.0,
    area=0.00785,
    density=7850.0,
    material_tag="steel_chain"
)
```

At this stage, you're working with **physical properties** (length, area, density, material).

### 2. **Geometry Generation Stage**

The package uses your `SegmentParameters` to create a discrete geometric model with mesh:

```julia
# Package internally creates geometry from your parameters
model = create_mooring_discrete_model(ph)
```

The length and point connections from `SegmentParameters` are used to generate the 1D mesh for each segment.

### 3. **FEM Setup Stage** (`MooringSegment`)

The package creates `MooringSegment` instances with **dedicated FE spaces** for each segment:

```julia
# Package internally creates MooringSegment from the discrete model
mooring_segment = MooringSegment(
    model,              # Discrete model with mesh
    "Segment_1",        # Tag from SegmentParameters
    pointA,             # MooringPoint (start)
    pointB,             # MooringPoint (stop)
    map_function,       # Geometric map (linear/quadratic)
    material,           # Material constitutive law
    density,            # Submerged density
    area                # Cross-sectional area
)

# Critical: Each segment gets its own FE space
U1, V1 = get_transient_FESpaces(segment1, order=1, dim=2)
U2, V2 = get_transient_FESpaces(segment2, order=1, dim=2)
# U1 and U2 are INDEPENDENT finite element spaces
```

At this stage, you're working with **FEM structures** (triangulations, FE spaces, residual functions).

### 4. **Why Each Segment Has Its Own FE Space**

This is a **key architectural feature** of Mooring.jl:

**Benefits:**
- **Flexible mesh refinement**: Each segment can have different element sizes
- **Independent DOFs**: Degrees of freedom are local to each segment
- **Easy assembly**: Segments are coupled only at shared points through boundary conditions
- **Modular design**: Add/remove segments without affecting others
- **Parallel potential**: Segments can be processed independently

**Example:**
```julia
# Three segments with different properties
segment1: order=1, 10 elements  → FE space V1 (11 nodes × dim DOFs)
segment2: order=2, 20 elements  → FE space V2 (41 nodes × dim DOFs)  
segment3: order=1, 5 elements   → FE space V3 (6 nodes × dim DOFs)

# Each has independent DOFs, coupled at endpoints
```

### 5. **What Each Structure Contains**

| Aspect | SegmentParameters | MooringSegment |
|--------|------------------|----------------|
| **Purpose** | User configuration | FEM computation |
| **Contains** | Length, area, density, material tag | Triangulation, FE space, material law, residual |
| **Mesh info** | None (just length) | Complete 1D triangulation |
| **FE space** | N/A | Dedicated TestFESpace & TrialFESpace |
| **Coupling** | Point IDs (start/stop) | MooringPoint instances with boundary conditions |
| **When used** | Parameter setup, YAML I/O | Assembly, solving, post-processing |
| **Typical user** | You (simulation setup) | Package internals (solver) |

### 6. **When You Interact With Each**

**Work with `SegmentParameters` when:**
- Reading/writing YAML configuration files
- Setting up simulation parameters (length, area, density)
- Defining material and drag tags
- Connecting segments to points
- Parametric studies (changing properties)

**Work with `MooringSegment` when:**
- Extending the package with custom residuals
- Implementing custom material models
- Developing advanced FEM features (higher-order elements)
- Debugging FEM assembly and solving
- Post-processing stress/strain fields

**Key Takeaway:** You primarily work with `SegmentParameters` for configuration. The package automatically handles:
- Conversion to `MooringSegment` with dedicated FE spaces
- Mesh generation for each segment
- Coupling between segments at shared points
- Assembly of the global FEM system

---

## Physical Map Functions

Mooring.jl provides functions to map from reference (1D) coordinates to physical (2D/3D) space:

### Linear Map

Creates a straight line between endpoints:
```julia
# Used automatically when segment length ≈ straight distance
# Good for taut cables with minimal sag
```

### Quadratic Map

Creates a catenary-like parabolic curve with correct arc length:
```julia
# Used automatically when segment length > straight distance
# Provides better initial guess for hanging cables
# Reduces Newton iterations in quasi-static analysis
```

The choice of map is handled automatically based on the geometry. See [MooringLine API](../API/Entities/MooringLine.md) for details.

---

## Best Practices

### Segment Length Guidelines

```julia
# Calculate straight-line distance
dx = x2 - x1
dy = y2 - y1  # (0 for 2D)
dz = z2 - z1
straight_distance = sqrt(dx^2 + dy^2 + dz^2)

# Recommended excess for catenary
if cable_type == "chain"
    excess = 0.05 to 0.15  # 5-15% excess
elseif cable_type == "synthetic"
    excess = 0.10 to 0.20  # 10-20% excess (more flexible)
end

segment_length = straight_distance * (1 + excess)
```

### Material Selection

| Application | Recommended Material | Why |
|-------------|---------------------|-----|
| Deep water (>1000m) | Polyester rope | Lighter weight, less bottom tension |
| Harsh environment | Chain | Abrasion resistance |
| Dynamic loads | Nylon rope | Energy absorption |
| Shallow water | Chain | Cost-effective, proven |

### Area and Density

Always verify your values:
```julia
# For circular cross-section
diameter = 0.120  # 120mm
area = π * (diameter/2)^2  # = 0.0113 m²

# Check submerged weight makes sense
weight_per_meter = (density - 1025.0) * area * 9.81
println("Submerged weight: $weight_per_meter N/m")
# Should be positive for steel, negative for buoyant synthetic ropes
```

---

## Common Issues

### Problem: Segment length too short

**Symptoms**: Solver fails to converge, unphysical tension
```julia
# Bad: Length less than straight distance
straight_dist = 100.0
segment_length = 90.0  # Too short!
```

**Solution**: Ensure length ≥ straight distance + sag allowance
```julia
# Good: Add appropriate excess
straight_dist = 100.0
segment_length = 110.0  # 10% excess for catenary
```

### Problem: Inconsistent units

**Symptoms**: Extremely large/small values, numerical issues
```julia
# Bad: Mixed units
area = 120.0       # Thought mm², actually m²!
density = 7.85     # Thought kg/m³, actually g/cm³!
```

**Solution**: Always use SI units
```julia
# Good: Consistent SI units
area = 0.0113      # m² (from 120mm diameter)
density = 7850.0   # kg/m³
```

### Problem: Wrong material tag

**Symptoms**: Error "material tag not found"
```julia
# Define material before referencing it
ph.materials["chain"] = PH.MaterialParameters(...)
ph.segments[1] = PH.SegmentParameters(material_tag="chain", ...)
```

---

## Next Steps

After defining segments:
1. Assemble segments into [Complete Lines](lines.md)
2. Configure [Material Properties](../API/Physics/Materials.md) in detail
3. Run [Quasi-Static Analysis](lines.md#solving-quasi-static-equilibrium)

## See Also

- [API: MooringSegment](../API/Entities/MooringSegment.md)
- [API: Materials](../API/Physics/Materials.md)
- [API: Drag Models](../API/Physics/Drag.md)
- [API: SeaBed Interaction](../API/Physics/SeaBed.md)