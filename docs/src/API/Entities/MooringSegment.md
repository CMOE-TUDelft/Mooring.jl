# MooringSegment

The `MooringSegment` module represents individual segments of a mooring line,
containing the physical properties and finite element discretization.

## Types

```@docs
Mooring.MooringSegments.MooringSegment
```

## Accessor Functions

```@docs
Mooring.MooringSegments.get_points
Mooring.MooringSegments.get_map
Mooring.MooringSegments.get_material
Mooring.MooringSegments.get_density
Mooring.MooringSegments.get_area
Mooring.MooringSegments.get_triangulation
Mooring.MooringSegments.get_measures
```

## Analysis Functions

```@docs
Mooring.MooringSegments.get_quasi_static_residual
```

Additional functions:
- `get_boundary_triangulation`: Returns the boundary triangulation
- `get_FE_spaces`: Returns finite element spaces for the segment

## Related

See [Material](../Physics/Materials.md) for material properties.

