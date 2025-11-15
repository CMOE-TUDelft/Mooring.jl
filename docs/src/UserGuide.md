# User Guide

This comprehensive user guide explains how to use [**Mooring.jl**](https://github.com/CMOE-TUDelft/Mooring.jl), a Julia package for high-fidelity modeling and simulation of mooring systems using the Finite Element Method.

## Overview

**Mooring.jl** provides a complete framework for simulating mooring lines with:
- **Finite Element Method (FEM)** with arbitrary order interpolation
- **Dynamic finite strain theory** for large deformation analysis
- **Tangential differential calculus** for geometric consistency
- **Linear and nonlinear material models** including viscoelasticity
- **Efficient time integration** for dynamic analysis
- Support for both **chain/steel cables** and **synthetic ropes**

## Structure of This Guide

This guide follows a **bottom-up approach**, building your understanding from fundamental components to complete simulations:

### 1. Setup and Installation
Start by installing Mooring.jl and understanding the configuration system.

### 2. Building Blocks
Learn about the core components:
- **Points**: Anchors, fairleads, and connection points with prescribed motions
- **Segments**: Cable sections with material properties and geometry
- **Lines**: Complete mooring lines connecting multiple segments

### 3. Physical Models
Configure the physics of your simulation:
- **Materials**: Elastic and viscoelastic constitutive models
- **Drag**: Hydrodynamic forces on the cables
- **Seabed**: Contact and interaction with the seafloor
- **Environmental Conditions**: Wave loading and current effects

### 4. Analysis and Results
Run simulations and analyze results:
- **Quasi-static Analysis**: Find equilibrium configurations
- **Dynamic Analysis**: Time-domain simulations
- **Visualization**: Plot and export results

## Learning Path

```@contents
Pages = [
    "userguide/installation.md",
    "userguide/points.md",
    "userguide/segments.md",
    "userguide/lines.md",
]
Depth = 1
```

## Quick Start Example

Here's a minimal example to get you started with a simple hanging cable:

```julia
using Mooring
import Mooring.ParameterHandlers as PH
import Mooring.MooringLines as ML

# Create parameter handler
ph = PH.ParameterHandler()

# Define material (steel chain)
ph.materials["chain"] = PH.MaterialParameters(
    tag="chain", 
    E=2.1e11,    # Young's modulus [Pa]
    μ=8.1e10     # Shear modulus [Pa]
)

# Define motion types
ph.motions["fixed"] = PH.MotionParameters(tag="fixed")

# Define points (anchor and fairlead)
ph.points[1] = PH.PointParameters(
    id=1, 
    coords=[0.0, 0.0, -100.0],  # Anchor at seabed
    motion_tag="fixed", 
    mesh_size=10.0
)
ph.points[2] = PH.PointParameters(
    id=2, 
    coords=[0.0, 0.0, -20.0],   # Fairlead at vessel
    motion_tag="fixed", 
    mesh_size=10.0
)

# Define segment
ph.segments[1] = PH.SegmentParameters(
    id=1,
    start_point=1,
    stop_point=2,
    material_tag="chain",
    area=0.01,           # Cross-sectional area [m²]
    density=7850.0,      # Density [kg/m³]
    length=120.0,        # Unstretched length [m]
    drag_tag="no_drag",
    seabed_tag="default_seabed"
)

# Define line
ph.lines[1] = PH.LineParameters(
    id=1,
    points=[1, 2],
    segments=[1]
)

# Solve quasi-static equilibrium
sol = ML.solve_quasistatic(ph)
```

## Key Concepts

### Parameter Handler
The `ParameterHandler` is the central configuration object that contains all parameters for your simulation:
- Material properties
- Point locations and motion types
- Segment geometry and properties
- Line topology

### Coordinate System
Mooring.jl uses a standard right-handed Cartesian coordinate system:
- **x, y**: Horizontal plane
- **z**: Vertical (positive upward, negative downward)
- Typical convention: z=0 at water surface, seabed at negative z

### Units
All quantities must be provided in **SI units**:
- Length: meters [m]
- Force: Newtons [N]
- Stress: Pascals [Pa]
- Density: kg/m³
- Time: seconds [s]

### Mesh Refinement
The `mesh_size` parameter at points controls the finite element mesh density:
- Smaller values → finer mesh → higher accuracy, longer computation
- Larger values → coarser mesh → lower accuracy, faster computation
- Typical range: 1-20 meters depending on problem scale

## Next Steps

1. **[Installation](userguide/installation.md)**: Set up Mooring.jl on your system
2. **[Points](userguide/points.md)**: Learn about point entities and motion types
3. **[Segments](userguide/segments.md)**: Configure cable segments with materials
4. **[Lines](userguide/lines.md)**: Assemble complete mooring lines

For API documentation, see the [API Reference](../API/API.md).

## Getting Help

If you encounter issues or have questions:
- Check the [API documentation](../API/API.md) for detailed function references
- Review the [examples](https://github.com/CMOE-TUDelft/Mooring.jl/tree/main/scripts)
- Open an [issue](https://github.com/CMOE-TUDelft/Mooring.jl/issues) on GitHub