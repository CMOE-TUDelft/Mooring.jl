# Defining Lines

Mooring lines are **assemblies of mooring segments and points** that define the geometry and behavior of a complete line in the mooring system.  
Each line connects (multiple) anchors to (multiple) fairleads, or intermediate points, through one or more segments.

## Table of Contents

- [Two Types of Line Structures](#two-types-of-line-structures)
- [Line Properties (LineParameters)](#line-properties-lineparameters)
- [Defining Lines in YAML](#defining-lines-in-yaml)
- [Defining Lines in Julia](#defining-lines-in-julia)
- [Workflow: From LineParameters to MooringLine](#workflow-from-lineparameters-to-mooringline)
- [Multi-Segment Lines](#multi-segment-lines)
- [Independent Line Solves](#independent-line-solves)
- [Best Practices](#best-practices)
- [See Also](#see-also)

## Two Types of Line Structures

**Mooring.jl** uses two distinct structures for handling lines at different levels:

### 1. `LineParameters` (High-Level: Parameter Handling)

`LineParameters` is used for **defining and configuring** lines in YAML files or through the parameter handler. This struct stores the topology of a line: which points and segments make up the line. You interact with this struct when setting up your mooring system configuration.

```julia
using Mooring.ParameterHandlers as PH

# Create line parameters for configuration
line_params = PH.LineParameters(
    id=1,
    points=[1, 2, 3],      # IDs of points that belong to this line
    segments=[1, 2]        # IDs of segments that form this line
)
```

### 2. `MooringLine` (Low-Level: FEM Operations)

`MooringLine` is used internally for **finite element analysis** operations. This struct contains a collection of `MooringSegment` instances and, critically, **creates a multifield FE space** by combining the individual FE spaces from each segment. The package automatically creates `MooringLine` instances from your `LineParameters` when setting up the simulation.

```julia
using Mooring.MooringLines

# MooringLine is created internally from the parameter handler
# You typically don't create these directly
mooring_line = MooringLine(segments_dict)

# Key feature: MooringLine creates a MultiFieldFESpace from segment FE spaces
X, Y = get_transient_FE_spaces(mooring_line)
# X is a TransientMultiFieldFESpace containing all segment trial spaces
# Y is a MultiFieldFESpace containing all segment test spaces
```

**Key Distinction:**
- **`LineParameters`**: What you define (point IDs, segment IDs) → *User configuration level*
- **`MooringLine`**: What the solver uses (multifield FE spaces, coupled residuals) → *FEM computation level*

**Critical Features:**
1. **Multifield FE Space**: `MooringLine` combines individual segment FE spaces into a single multifield space
2. **Independent Solve**: Each `MooringLine` results in an **independent solve statement** (separate nonlinear system)
3. **Line-Level Coupling**: Segments within a line are strongly coupled through the multifield formulation
4. **Inter-Line Independence**: Different lines do not interact (can be solved in parallel)

## Line Properties (LineParameters)

When defining lines through the parameter handler, each line uses the following parameters:

| Parameter | Type | Description |
|-----------|------|-------------|
| `id` | `Int` | Unique integer identifier for the line |
| `tag` | `String` | Human-readable name for the line (auto-generated if not provided) |
| `points` | `Vector{Int}` | Ordered list of point IDs that belong to this line |
| `segments` | `Vector{Int}` | List of segment IDs that form the line |

### Points Vector

The `points` vector defines all points (anchors, fairleads, intermediate) that belong to the line:
- Must be ordered from one end to the other
- Must include all endpoints of all segments in the line
- Typically: `[anchor_id, intermediate_1, intermediate_2, ..., fairlead_id]`

### Segments Vector

The `segments` vector lists all segments that compose the line:
- Order doesn't matter for the solver (topology is determined from segment start/stop points)
- Each segment must connect two points from the `points` vector
- Segments are automatically coupled where they share common points

---

## Defining Lines in YAML

### Single Segment Line

```yaml
# Define points first
points:
  - id: 1
    coords: [0.0, 0.0, -50.0]
    motion_tag: "fixed"
    mesh_size: 10.0
  - id: 2
    coords: [10.0, 0.0, 0.0]
    motion_tag: "fairlead_motion"
    mesh_size: 5.0

# Define segments
segments:
  - id: 1
    start_point: 1
    stop_point: 2
    length: 55.0
    area: 0.01
    density: 7850.0
    material_tag: "steel"

# Define line connecting them
lines:
  - id: 1
    points: [1, 2]
    segments: [1]
```

In this example:
- A single line is created with two points (anchor at seabed, fairlead at floater)
- Connected by one segment
- Results in **one solve statement** with a single-field FE space (only one segment)

### Multi-Segment Line (Hybrid)

```yaml
points:
  - id: 1
    coords: [0.0, 0.0, -100.0]
    motion_tag: "fixed"
  - id: 2
    coords: [25.0, 0.0, -50.0]
    motion_tag: "fixed"
  - id: 3
    coords: [50.0, 0.0, 0.0]
    motion_tag: "fairlead_motion"

segments:
  - id: 1
    start_point: 1
    stop_point: 2
    length: 55.0
    material_tag: "chain"
    area: 0.007
    density: 7700.0
  - id: 2
    start_point: 2
    stop_point: 3
    length: 60.0
    material_tag: "polyester"
    area: 0.004
    density: 1400.0

lines:
  - id: 1
    points: [1, 2, 3]
    segments: [1, 2]
```

This defines a hybrid line:
- Lower segment: chain (anchor → midpoint)
- Upper segment: polyester rope (midpoint → fairlead)
- Results in **one solve statement** with a **multifield FE space** (two segment spaces combined)

---

## Defining Lines in Julia

### Basic Single-Segment Line

```julia
using Mooring
import Mooring.ParameterHandlers as PH

# Initialize parameter handler
ph = PH.ParameterHandler()

# Define points
ph.points[1] = PH.PointParameters(
    id=1, 
    coords=[0.0, 0.0, -50.0], 
    motion_tag="fixed"
)
ph.points[2] = PH.PointParameters(
    id=2, 
    coords=[10.0, 0.0, 0.0], 
    motion_tag="fairlead_motion"
)

# Define material
ph.materials["steel"] = PH.MaterialParameters(
    E=2.1e11,  # Young's modulus [Pa]
    ν=0.3      # Poisson's ratio
)

# Define segment
ph.segments[1] = PH.SegmentParameters(
    id=1, 
    start_point=1, 
    stop_point=2,
    material_tag="steel", 
    area=0.01,
    density=7850.0, 
    length=55.0
)

# Define line
ph.lines[1] = PH.LineParameters(
    id=1, 
    points=[1, 2], 
    segments=[1]
)
```

### Multi-Segment Hybrid Line

```julia
using Mooring
import Mooring.ParameterHandlers as PH

ph = PH.ParameterHandler()

# Define points (anchor, midpoint, fairlead)
ph.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0, -100.0], motion_tag="fixed")
ph.points[2] = PH.PointParameters(id=2, coords=[25.0, 0.0, -50.0], motion_tag="fixed")
ph.points[3] = PH.PointParameters(id=3, coords=[50.0, 0.0, 0.0], motion_tag="fairlead_motion")

# Define materials
ph.materials["chain"] = PH.MaterialParameters(E=2.1e11, ν=0.3)
ph.materials["polyester"] = PH.MaterialParameters(E=5.0e9, ν=0.3)

# Define segments
ph.segments[1] = PH.SegmentParameters(
    id=1, 
    start_point=1, 
    stop_point=2,
    length=55.0,
    material_tag="chain",
    area=0.007,
    density=7700.0
)

ph.segments[2] = PH.SegmentParameters(
    id=2, 
    start_point=2, 
    stop_point=3,
    length=60.0,
    material_tag="polyester",
    area=0.004,
    density=1400.0
)

# Define line with two segments
ph.lines[1] = PH.LineParameters(
    id=1, 
    points=[1, 2, 3],    # All three points
    segments=[1, 2]       # Both segments
)
```

---

## Workflow: From LineParameters to MooringLine

Understanding how `LineParameters` transforms into `MooringLine` reveals the multifield FE space architecture:

### 1. **User Configuration Stage** (`LineParameters`)

You define the line topology by specifying which points and segments belong together:

```julia
# Define configuration - just topology (IDs)
ph.lines[1] = PH.LineParameters(
    id=1,
    points=[1, 2, 3],
    segments=[1, 2]
)
```

At this stage, you're working with **topology** (which points/segments form a line).

### 2. **Setup Stage: Creating MooringLine**

The package creates a `MooringLine` by:
1. Creating a discrete model for the line
2. Instantiating `MooringPoint` objects from point IDs
3. Instantiating `MooringSegment` objects from segment IDs
4. Storing segments in a dictionary

```julia
# Package internally creates MooringLine from LineParameters
# This happens in setup_lines(ph)

# For each line:
model = generate_discrete_model(line_params, ph)

# Create MooringPoints
points = Dict{Int, MooringPoint}()
for p_id in line_params.points
    points[p_id] = MooringPoint(model, point_params.tag, motion)
end

# Create MooringSegments (each with its own FE space)
segments = Dict{Int, MooringSegment}()
for s_id in line_params.segments
    segment = MooringSegment(...)
    segments[s_id] = segment
end

# Create MooringLine
mooring_line = MooringLine(segments)
```

### 3. **FEM Setup Stage: Multifield FE Space Construction**

This is where the **critical transformation** happens. `MooringLine` creates a **multifield FE space**:

```julia
# Get transient FE spaces for the line
X, Y = get_transient_FE_spaces(mooring_line)

# Internally, this does:
# 1. Get individual FE spaces from each segment
test_spaces = []
trial_spaces = []
for segment in mooring_line.segments
    U_seg, V_seg = get_transient_FESpaces(segment)
    push!(trial_spaces, U_seg)  # Trial space for this segment
    push!(test_spaces, V_seg)   # Test space for this segment
end

# 2. Combine into multifield spaces
X = TransientMultiFieldFESpace(trial_spaces)  # All trial spaces
Y = MultiFieldFESpace(test_spaces)            # All test spaces
```

**Result**: 
- `X` and `Y` are **multifield spaces** containing all segment DOFs
- Segments are coupled through shared boundary points
- This creates a **single unified FE space for the entire line**

### 4. **Why Multifield FE Spaces?**

The multifield FE space architecture provides:

**Benefits:**
- **Modular assembly**: Each segment contributes its own residual
- **Flexible coupling**: Segments coupled only at shared points
- **Independent segment properties**: Different mesh, order, material per segment
- **Efficient assembly**: Leverage Gridap's multifield infrastructure
- **Clear structure**: Each field corresponds to one segment

**Example with 2 segments:**
```julia
# Segment 1: 10 elements, order 1 → 11 nodes × 2D = 22 DOFs
# Segment 2: 20 elements, order 1 → 21 nodes × 2D = 42 DOFs

# Multifield space:
# - Field 1: 22 DOFs (segment 1)
# - Field 2: 42 DOFs (segment 2)
# - Total: 64 DOFs before applying boundary conditions
# - Coupling: Enforced at shared point (end of seg1 = start of seg2)
```

### 5. **Independent Solve per Line**

Each `MooringLine` results in an **independent nonlinear solve**:

```julia
# For each line in the mooring system:
for (line_id, line) in mooring_lines
    
    # Get multifield FE spaces for THIS line
    X, Y = get_transient_FE_spaces(line)
    
    # Get reference configuration for THIS line
    Xₕ = get_reference_configuration(line, X(0.0))
    
    # Get residual for THIS line (sum of segment residuals)
    res = get_quasi_static_residual(line, Xₕ)
    
    # Solve THIS line independently
    op = FEOperator(res, X(0.0), Y)
    nls = NLSolver(...)
    uₕ = solve(nls, op)  # ← Independent solve for this line
    
    # Store solution for this line
    solutions[line_id] = uₕ
end
```

**Key Points:**
- **One solve per line**: Each line is solved independently
- **No inter-line coupling**: Lines don't interact in the solver
- **Parallel potential**: Different lines can be solved on different processors
- **Separate convergence**: Each line converges independently

### 6. **What Each Structure Contains**

| Aspect | LineParameters | MooringLine |
|--------|---------------|-------------|
| **Purpose** | User configuration | FEM computation |
| **Contains** | Point IDs, segment IDs | MooringSegment instances |
| **FE space** | N/A | MultiFieldFESpace (combines segment spaces) |
| **Topology** | Explicit (IDs) | Implicit (through segment connections) |
| **Solve** | N/A | One independent solve per line |
| **Coupling** | None (just IDs) | Segments coupled within line, lines independent |
| **When used** | Parameter setup, YAML I/O | Assembly, solving, post-processing |
| **Typical user** | You (simulation setup) | Package internals (solver) |

**Key Takeaway:** 
- You define lines using `LineParameters` (just topology)
- Package creates `MooringLine` with multifield FE spaces (one field per segment)
- Each line is solved independently (one solve statement per line)
- Segments within a line are strongly coupled; lines are independent

---

## Multi-Segment Lines

Multi-segment lines are common in practical mooring systems:

### Chain-Rope-Chain Configuration

```julia
# Ground chain + polyester rope + top chain
ph.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0, -500.0], motion_tag="fixed")
ph.points[2] = PH.PointParameters(id=2, coords=[200.0, 0.0, -480.0], motion_tag="fixed")
ph.points[3] = PH.PointParameters(id=3, coords=[800.0, 0.0, -100.0], motion_tag="fixed")
ph.points[4] = PH.PointParameters(id=4, coords=[900.0, 0.0, -20.0], motion_tag="vessel_motion")

# Bottom chain
ph.segments[1] = PH.SegmentParameters(
    id=1, start_point=1, stop_point=2, length=220.0,
    material_tag="chain", area=0.0113, density=7700.0
)

# Polyester rope (lighter, more compliant)
ph.segments[2] = PH.SegmentParameters(
    id=2, start_point=2, stop_point=3, length=850.0,
    material_tag="polyester", area=0.0314, density=1100.0
)

# Top chain
ph.segments[3] = PH.SegmentParameters(
    id=3, start_point=3, stop_point=4, length=120.0,
    material_tag="chain", area=0.0113, density=7700.0
)

# One line with three segments
ph.lines[1] = PH.LineParameters(
    id=1, 
    points=[1, 2, 3, 4],      # Four points
    segments=[1, 2, 3]         # Three segments
)

# This creates:
# - One MooringLine
# - One multifield FE space with 3 fields (one per segment)
# - One solve statement
```

### Advantages of Multi-Segment Lines

1. **Material variation**: Different properties along the line
2. **Cost optimization**: Heavy chain near seabed, lighter rope in water column
3. **Dynamic response**: Tuned compliance at different depths
4. **Mesh refinement**: Finer mesh in high-stress regions
5. **Physical accuracy**: Model actual mooring compositions

---

## Independent Line Solves

### Why Each Line is Solved Independently

**Physical reasoning:**
- Mooring lines typically don't touch each other
- Each line has independent boundary conditions
- Load on one line doesn't directly affect others (decoupled)

**Computational benefits:**
- **Smaller systems**: Each line is smaller than global system
- **Better convergence**: Localized nonlinearities easier to resolve
- **Parallel execution**: Lines can be solved simultaneously
- **Modular debugging**: Test individual lines in isolation
- **Flexible analysis**: Solve only selected lines

### Solve Structure

```julia
function solve_quasistatic(ph::ParameterHandler)
    # Setup all lines
    mooring_lines = setup_lines(ph)
    
    # Storage for solutions
    solutions = Dict{Int, FEFunction}()
    
    # Solve each line independently
    for (line_id, line) in mooring_lines
        @info "Solving line $line_id with $(length(line.segments)) segments"
        
        # Create multifield FE spaces for this line
        X, Y = get_transient_FE_spaces(line)
        
        # Reference configuration
        Xₕ = get_reference_configuration(line, X(0.0))
        
        # Assemble residual (sum of segment residuals)
        res = get_quasi_static_residual(line, Xₕ)
        
        # Define FE operator
        op = FEOperator(res, X(0.0), Y)
        
        # Solve nonlinear system for this line
        nls = NLSolver(
            BackslashSolver(), 
            iterations=200, 
            show_trace=true, 
            ftol=1e-8,
            method=:newton
        )
        uₕ = solve(nls, op)
        
        # Store solution
        solutions[line_id] = uₕ
    end
    
    return solutions
end
```

### Multi-Line System Example

```yaml
# Three independent mooring lines (spread mooring)
lines:
  - id: 1  # Line at 0° azimuth
    points: [1, 2, 3]
    segments: [1, 2]
    
  - id: 2  # Line at 120° azimuth
    points: [4, 5, 6]
    segments: [3, 4]
    
  - id: 3  # Line at 240° azimuth
    points: [7, 8, 9]
    segments: [5, 6]
```

**Solver behavior:**
```
Solving line 1 with 2 segments...
  Newton iteration 1: |res| = 1.2e+05
  Newton iteration 2: |res| = 3.4e+03
  ...
  Newton iteration 15: |res| = 8.7e-09  ✓ Converged

Solving line 2 with 2 segments...
  Newton iteration 1: |res| = 1.1e+05
  ...
  Newton iteration 12: |res| = 5.2e-09  ✓ Converged

Solving line 3 with 2 segments...
  ...
  ✓ Converged
```

Each line converges independently!

---

## Best Practices

### Line Definition

1. **Order points correctly**: List points in order along the line
   ```julia
   # Good: ordered from anchor to fairlead
   points=[1, 2, 3, 4]  # anchor → mid1 → mid2 → fairlead
   
   # Also OK: segments determine actual topology
   points=[1, 3, 2, 4]  # Order doesn't affect solver, but less clear
   ```

2. **Include all necessary points**: Every segment endpoint must be in the points list
   ```julia
   # Bad: missing intermediate point
   points=[1, 3]      # Where is point 2?
   segments=[1, 2]    # But segment 1 ends at point 2!
   
   # Good: all points included
   points=[1, 2, 3]
   segments=[1, 2]
   ```

3. **Use descriptive tags** (optional but recommended):
   ```julia
   ph.lines[1] = PH.LineParameters(
       id=1,
       tag="starboard_forward_line",  # Clear identification
       points=[1, 2, 3],
       segments=[1, 2]
   )
   ```

### Multi-Line Systems

1. **Ensure lines are physically separated**: Lines should not share segments
   ```julia
   # Good: independent lines
   ph.lines[1] = PH.LineParameters(id=1, segments=[1, 2])
   ph.lines[2] = PH.LineParameters(id=2, segments=[3, 4])
   
   # Bad: lines share segment 2
   ph.lines[1] = PH.LineParameters(id=1, segments=[1, 2])
   ph.lines[2] = PH.LineParameters(id=2, segments=[2, 3])  # Segment 2 used twice!
   ```

2. **Lines can share fairlead points** (same fairlead, different lines):
   ```julia
   # OK: different lines can connect to the same fairlead
   ph.points[10] = PH.PointParameters(id=10, coords=[0.0, 0.0, 0.0], motion_tag="vessel")
   ph.lines[1] = PH.LineParameters(id=1, points=[1, 2, 10], segments=[1, 2])
   ph.lines[2] = PH.LineParameters(id=2, points=[3, 4, 10], segments=[3, 4])
   # Both lines end at point 10, but solve independently
   ```

### Performance Considerations

1. **Balance segment count**: More segments = more fields = larger multifield space
   ```julia
   # Trade-off: accuracy vs. computational cost
   # More segments: better material variation representation, larger system
   # Fewer segments: faster solve, less flexibility
   ```

2. **Symmetry**: Identical lines can reuse geometry/materials
   ```julia
   # Define material once, use for multiple segments/lines
   ph.materials["chain"] = PH.MaterialParameters(...)
   ph.segments[1] = PH.SegmentParameters(..., material_tag="chain")
   ph.segments[2] = PH.SegmentParameters(..., material_tag="chain")
   ```

---

## See Also

- [API: MooringLine](../API/Entities/MooringLine.md)
- [User Guide: Points](points.md)
- [User Guide: Segments](segments.md)
- [API: MultiField FE Spaces](https://gridap.github.io/Gridap.jl/stable/MultiField/)