# MooringLine

The `MooringLine` module provides functionality for managing complete mooring lines,
which consist of one or more segments connecting multiple points.

## Types

```@docs
Mooring.MooringLines.MooringLine
```

## Functions

```@docs
Mooring.MooringLines.get_segments
Mooring.MooringLines.setup_lines
Mooring.MooringLines.get_transient_FE_spaces
Mooring.MooringLines.get_reference_configuration
Mooring.MooringLines.get_quasi_static_residual
Mooring.MooringLines.solve_quasistatic
```

## Physical Map Functions

These functions create mappings from reference coordinates to physical coordinates:

- `get_physical_linear_map`: Creates a linear interpolation between endpoints
- `get_physical_quadratic_map`: Creates a catenary-like parabolic curve with correct arc length

See [`MooringSegment.get_quasi_static_residual`](@ref) for related segment-level functionality.

