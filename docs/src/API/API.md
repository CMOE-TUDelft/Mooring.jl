```@meta
CurrentModule = Mooring
```

# API Documentation

This section provides detailed API documentation for all modules in Mooring.jl.

## Modules Overview

### Entities
Core data structures representing mooring system components:
- [MooringLine](Entities/MooringLine.md) - Line entities connecting points
- [MooringPoint](Entities/MooringPoint.md) - Point entities (anchors, fairleads, etc.)
- [MooringSegment](Entities/MooringSegment.md) - Segment entities within lines

### Geometry
Geometric mesh generation and discrete models:
- [MooringDiscreteModel](Geometry/MooringDiscreteModel.md) - Discrete model generation

### IO
Input/output and parameter handling:
- [ParameterHandler](IO/ParameterHandler.md) - Configuration and parameter management

### Physics
Physical models and material properties:
- [Materials](Physics/Materials.md) - Material constitutive models
- [Drag](Physics/Drag.md) - Hydrodynamic drag models
- [EnvironmentalConditions](Physics/EnvironmentalConditions.md) - Environmental loading
- [PointMotion](Physics/PointMotion.md) - Prescribed motion at points
- [SeaBed](Physics/SeaBed.md) - Seabed contact and interaction
- [TangentialDiffCalculus](Physics/TangentialDiffCalculus.md) - Differential geometry operations