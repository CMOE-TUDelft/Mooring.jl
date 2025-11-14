# Materials

The `Materials` module defines material constitutive models for mooring cables.

## Types

```@docs
Mooring.Materials.Material
```

## Functions

Accessor functions for material properties:
- `get_young_modulus`: Returns Young's modulus (E)
- `get_shear_modulus`: Returns shear modulus (μ)  
- `get_lame_lambda`: Returns Lamé's first parameter (λ)

These properties are used in the constitutive relations for cable mechanics.

