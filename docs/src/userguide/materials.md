# Defining Materials

Materials define the constitutive behavior of mooring segments, describing how cables respond to tension and deformation. **Mooring.jl** provides material models for both elastic and viscoelastic* behavior (*under development).

## Table of Contents

- [Two Types of Material Structures](#two-types-of-material-structures)
- [Material Properties (MaterialParameters)](#material-properties-materialparameters)
- [Available Material Models](#available-material-models)
- [Defining Materials in YAML](#defining-materials-in-yaml)
- [Defining Materials in Julia](#defining-materials-in-julia)
- [Workflow: From MaterialParameters to Material](#workflow-from-materialparameters-to-material)
- [Stress Measures](#stress-measures)
- [Material Selection Guide](#material-selection-guide)
- [Best Practices](#best-practices)
- [See Also](#see-also)

## Two Types of Material Structures

**Mooring.jl** uses two distinct structures for handling materials at different levels:

### 1. `MaterialParameters` (High-Level: Parameter Handling)

`MaterialParameters` is used for **defining and configuring** material properties in YAML files or through the parameter handler. This struct stores all the user-specified properties needed to create a material model. You interact with this struct when setting up your simulation configuration.

```julia
using Mooring.ParameterHandlers as PH

# Create material parameters for configuration
mat_params = PH.MaterialParameters(
    tag="steel_chain",
    type="LinearElastic",
    E=2.1e11,      # Young's modulus [Pa]
    μ=8.08e10      # Shear modulus [Pa]
)
```

### 2. `Material` (Low-Level: FEM Operations)

`Material` is an abstract type with concrete implementations (`LinearElastic`, `Scharpery`) used internally for **finite element analysis** operations. These structs contain the constitutive model and stress calculation functions needed for FEM computations. The package automatically creates `Material` instances from your `MaterialParameters`.

```julia
using Mooring.Materials

# Material is created internally from MaterialParameters
# You typically don't create these directly
material = LinearElastic(E=2.1e11, μ=8.08e10)

# Used to compute stresses during FEM solve
S_stress = S(material, Etang)  # Second Piola-Kirchhoff stress
```

**Key Distinction:**
- **`MaterialParameters`**: What you define (properties, model type) → *User configuration level*
- **`Material`**: What the solver uses (constitutive laws, stress functions) → *FEM computation level*

**Critical Feature:** The `Material` structs implement stress calculation functions that are called during FEM assembly to compute internal forces in the mooring segments.

## Material Properties (MaterialParameters)

When defining materials through the parameter handler, each material uses the following base parameters:

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `tag` | `String` | Yes | Unique identifier for the material |
| `type` | `String` | Yes | Material model type: "LinearElastic" or "Scharpery" |

Additional parameters depend on the material type selected.

---

## Available Material Models

### 1. Linear Elastic Material

**Type**: `"LinearElastic"`

Models linear elastic behavior following Hooke's law. Suitable for steel chains, wire ropes, and materials with negligible time-dependent behavior.

**Parameters:**

| Parameter | Type | Default | Units | Description |
|-----------|------|---------|-------|-------------|
| `E` | `Real` | 1.0 | Pa | Young's modulus (stiffness) |
| `μ` | `Real` | 0.5*E | Pa | Second Lamé constant (shear modulus) |

**Constitutive Law:**

The second Piola-Kirchhoff stress is computed as:

$$\mathbf{S} = 2\mu \mathbf{E}_{\text{tang}}$$

where $\mathbf{E}_{\text{tang}}$ is the tangential Lagrangian Green strain.

**Typical Values:**

| Material | Young's Modulus E [GPa] | Shear Modulus μ [GPa] |
|----------|-------------------------|------------------------|
| Steel (chain/wire) | 200-210 | 77-81 |
| Polyester rope | 5-10 | 2-4 |
| Nylon rope | 2-5 | 0.8-2 |
| HMPE (high-modulus polyethylene) | 80-120 | 30-45 |

**Example:**

```julia
# Steel chain
ph.materials["steel"] = PH.MaterialParameters(
    tag="steel",
    type="LinearElastic",
    E=2.1e11,      # 210 GPa
    μ=8.08e10      # 80.8 GPa
)

# Polyester rope
ph.materials["polyester"] = PH.MaterialParameters(
    tag="polyester",
    type="LinearElastic",
    E=5.0e9,       # 5 GPa
    μ=1.92e9       # 1.92 GPa
)
```

---

### 2. Scharpery Viscoelastic Material (under development)

**Type**: `"Scharpery"`

Models nonlinear viscoelastic behavior using the Scharpery model. Suitable for synthetic fiber ropes (polyester, nylon) that exhibit time-dependent creep and stress relaxation.

**Parameters:**

| Parameter | Type | Default | Units | Description |
|-----------|------|---------|-------|-------------|
| `D0` | `Real` | 1.0 | Pa⁻¹ | Linear elastic compliance (instantaneous) |
| `N` | `Int` | 1 | - | Number of relaxation times |
| `Dn` | `Vector{Real}` | [1.0] | Pa⁻¹ | Compliance for each relaxation time |
| `λn` | `Vector{Real}` | [1.0] | s⁻¹ | Inverse of relaxation times |
| `g0` | `String` | "1.0" | - | Nonlinear coefficient for instantaneous compliance (function of σ) |
| `g1` | `String` | "1.0" | - | Nonlinear coefficient for transient compliance (function of σ) |
| `g2` | `String` | "1.0" | - | Nonlinear coefficient for stress rate-dependent compliance (function of σ) |

**Constitutive Law:**

The Scharpery model relates strain to stress history through:

$$\varepsilon(t) = g_0 D_0 \sigma(t) + g_1 \sum_{n=1}^{N} D_n \int_0^t e^{-\lambda_n (t-\tau)} g_2(\sigma(\tau)) \frac{d\sigma}{d\tau} d\tau$$

where:
- $D_0$: Instantaneous elastic compliance
- $D_n$: Transient compliance for each mode
- $\lambda_n$: Relaxation rate for each mode
- $g_0, g_1, g_2$: Nonlinear stress-dependent functions

**Example:**

```julia
# Polyester rope with single relaxation time
ph.materials["polyester_viscoelastic"] = PH.MaterialParameters(
    tag="polyester_viscoelastic",
    type="Scharpery",
    D0=2.0e-10,           # Instantaneous compliance [Pa⁻¹]
    N=1,                  # Single relaxation mode
    Dn=[1.0e-10],         # Transient compliance [Pa⁻¹]
    λn=[1.0e-3],          # Relaxation rate [s⁻¹]
    g0="1.0",             # Linear instantaneous response
    g1="1.0",             # Linear transient response
    g2="1.0"              # Linear stress-rate response
)

# Advanced: Multiple relaxation times
ph.materials["nylon_multimode"] = PH.MaterialParameters(
    tag="nylon_multimode",
    type="Scharpery",
    D0=5.0e-10,
    N=3,                  # Three relaxation modes
    Dn=[2.0e-10, 1.5e-10, 1.0e-10],
    λn=[1.0e-2, 1.0e-3, 1.0e-4],  # Fast, medium, slow relaxation
    g0="1.0 + 0.01*σ",    # Slightly nonlinear instantaneous
    g1="1.0 + 0.02*σ",    # More nonlinear transient
    g2="1.0"
)
```

---

## Defining Materials in YAML

### Simple Linear Elastic Material

```yaml
materials:
  - tag: steel_chain
    type: LinearElastic
    E: 2.1e11        # Young's modulus [Pa]
    μ: 8.08e10       # Shear modulus [Pa]
```

### Multiple Materials

```yaml
materials:
  # Steel chain for lower segment
  - tag: chain
    type: LinearElastic
    E: 2.1e11
    μ: 8.08e10
  
  # Polyester rope for upper segment
  - tag: polyester
    type: LinearElastic
    E: 5.0e9
    μ: 1.92e9
  
  # HMPE rope (high stiffness synthetic)
  - tag: hmpe
    type: LinearElastic
    E: 1.0e11
    μ: 3.85e10
```

### Viscoelastic Material

```yaml
materials:
  - tag: polyester_viscoelastic
    type: Scharpery
    D0: 2.0e-10       # Instantaneous compliance [Pa⁻¹]
    N: 2              # Two relaxation modes
    Dn: [1.5e-10, 1.0e-10]
    λn: [1.0e-2, 1.0e-3]
    g0: "1.0"
    g1: "1.0"
    g2: "1.0"
```

Load in Julia:
```julia
using Mooring
import Mooring.ParameterHandlers as PH

ph = PH.read_parameters("mooring_config.yaml")
```

---

## Defining Materials in Julia

### Basic Linear Elastic

```julia
using Mooring
import Mooring.ParameterHandlers as PH

# Initialize parameter handler
ph = PH.ParameterHandler()

# Define steel chain material
ph.materials["steel"] = PH.MaterialParameters(
    tag="steel",
    type="LinearElastic",
    E=2.1e11,      # 210 GPa
    μ=8.08e10      # 80.8 GPa (computed from E and Poisson's ratio ν=0.3)
)

# Define polyester rope material
ph.materials["polyester"] = PH.MaterialParameters(
    tag="polyester",
    type="LinearElastic",
    E=5.0e9,       # 5 GPa
    μ=1.92e9       # 1.92 GPa
)
```

### Computing Shear Modulus from Poisson's Ratio

```julia
# Material properties
E = 2.1e11         # Young's modulus [Pa]
ν = 0.3            # Poisson's ratio [-]

# Compute shear modulus
μ = E / (2 * (1 + ν))

ph.materials["steel"] = PH.MaterialParameters(
    tag="steel",
    type="LinearElastic",
    E=E,
    μ=μ
)
```

### Viscoelastic Material

```julia
# Polyester rope with time-dependent behavior
ph.materials["polyester_visco"] = PH.MaterialParameters(
    tag="polyester_visco",
    type="Scharpery",
    D0=2.0e-10,           # Instantaneous compliance
    N=1,                  # Single relaxation mode
    Dn=[1.0e-10],         # Transient compliance
    λn=[1.0e-3],          # Relaxation rate (1000 s characteristic time)
    g0="1.0",             # Linear behavior
    g1="1.0",
    g2="1.0"
)
```

---

## Workflow: From MaterialParameters to Material

Understanding how `MaterialParameters` transforms into `Material` helps clarify the constitutive modeling:

### 1. **User Configuration Stage** (`MaterialParameters`)

You define material properties in YAML or Julia:

```julia
# Define configuration
ph.materials["steel"] = PH.MaterialParameters(
    tag="steel",
    type="LinearElastic",
    E=2.1e11,
    μ=8.08e10
)
```

At this stage, you're working with **material properties** (Young's modulus, compliance, etc.).

### 2. **Segment Creation Stage**

When creating segments, you reference materials by their tag:

```julia
ph.segments[1] = PH.SegmentParameters(
    id=1,
    start_point=1,
    stop_point=2,
    material_tag="steel",  # Reference to material
    ...
)
```

### 3. **Material Instantiation Stage**

The package creates `Material` instances when building `MooringSegment` objects:

```julia
# Package internally creates Material from MaterialParameters
# This happens in setup_lines(ph)

mat_params = ph.materials["steel"]
material = Material(mat_params)

# For LinearElastic, this creates:
# material = LinearElastic(E=2.1e11, μ=8.08e10)
```

### 4. **FEM Assembly Stage: Stress Calculation**

The `Material` instances are used to compute stresses during FEM assembly:

```julia
# During residual assembly for a segment

# given tangential strain tensor is computed (see TangentialDiffCalculus.jl), compute second Piola-Kirchhoff stress using material law
S_stress = S(material, Etang)
# For LinearElastic: S = 2*μ*Etang

# Then, compute first Piola-Kirchhoff stress
K_stress = K(FΓ, S_stress)
# K = FΓ ⋅ S

# Finally, compute the residual
res = ∫((∇(v)' ⋅ Q') ⊙ K_stress * Jabs)dΩ
```

### 5. **What Each Structure Contains**

| Aspect | MaterialParameters | Material (LinearElastic/Scharpery) |
|--------|-------------------|-------------------------------------|
| **Purpose** | User configuration | FEM computation |
| **Contains** | Material properties (E, μ, D0, etc.) | Constitutive law implementation |
| **Functions** | None (just data) | `S(material, Etang)` - stress calculation |
| **When used** | Parameter setup, YAML I/O | Assembly, stress computation |
| **Typical user** | You (simulation setup) | Package internals (solver) |
| **Mutability** | Can be modified before setup | Fixed once created |

### 6. **When You Interact With Each**

**Work with `MaterialParameters` when:**
- Reading/writing YAML configuration files
- Defining material properties for segments
- Conducting parametric studies (varying E, μ)
- Selecting material models (LinearElastic vs Scharpery)

**Work with `Material` when:**
- Implementing custom constitutive models
- Extending the package with new material types
- Developing stress calculation methods
- Debugging constitutive behavior

### 7. **Example: Complete Workflow**

```julia
using Mooring
import Mooring.ParameterHandlers as PH

# Step 1: Define material parameters (MaterialParameters)
ph = PH.ParameterHandler()
ph.materials["steel"] = PH.MaterialParameters(
    tag="steel",
    type="LinearElastic",
    E=2.1e11,
    μ=8.08e10
)

# Step 2: Reference material in segment definition
ph.segments[1] = PH.SegmentParameters(
    id=1,
    start_point=1,
    stop_point=2,
    material_tag="steel",  # Links to material
    area=0.00785,
    density=7850.0,
    length=100.0
)

# Step 3: Create mooring lines (internally creates Material instances)
mooring_lines = setup_lines(ph)

# Internally, for each segment:
# mat_params = ph.materials["steel"]
# material = Material(mat_params)  # Creates LinearElastic(E=2.1e11, μ=8.08e10)
# segment = MooringSegment(..., material, ...)

# Step 4: Solve (Material used in stress calculations during assembly)
solutions = solve_quasistatic(ph)
```

**Key Takeaway:** You define materials using `MaterialParameters` (properties). The package creates `Material` instances (constitutive laws) for stress calculations during FEM assembly.

---

## Stress Measures

**Mooring.jl** computes several stress measures internally:

### Second Piola-Kirchhoff Stress ($\mathbf{S}$)

**Definition:**
$$\mathbf{S} = \frac{\partial W}{\partial \mathbf{E}_{\text{tang}}}$$

where $W$ is the strain energy density and $\mathbf{E}_{\text{tang}}$ is the tangential Lagrangian Green strain.

**For LinearElastic:**
$$\mathbf{S} = 2\mu \mathbf{E}_{\text{tang}}$$

**Used for:** Constitutive law implementation (energy-conjugate to strain)

### First Piola-Kirchhoff Stress ($\mathbf{K}$)

**Definition:**
$$\mathbf{K} = \mathbf{F}_\Gamma \cdot \mathbf{S}$$

where $\mathbf{F}_\Gamma$ is the deformation gradient along the line.

**Used for:** Weak form assembly in FEM residual

### Cauchy Stress ($\boldsymbol{\sigma}$)

**Definition:**
$$\boldsymbol{\sigma} = \frac{1}{\Lambda} \mathbf{K} \cdot \mathbf{F}_\Gamma^T$$

where $\Lambda$ is the stretch of the line.

**Used for:** Physical interpretation (true stress in deformed configuration), post-processing

---

## Material Selection Guide (not rigurously tested)

### By Application

| Application | Recommended Material | Type | Why |
|-------------|---------------------|------|-----|
| Permanent mooring (deep water) | Polyester rope | LinearElastic or Scharpery | Lightweight, cost-effective, good fatigue |
| Permanent mooring (harsh environment) | Chain | LinearElastic | Abrasion resistance, proven reliability |
| Taut-leg mooring | HMPE or polyester | LinearElastic | High stiffness-to-weight ratio |
| Compliant mooring | Nylon or polyester | Scharpery | Energy absorption, compliance |
| Temporary mooring | Wire rope | LinearElastic | High strength, low cost |
| Ultra-deep water (>2000m) | Synthetic fiber | LinearElastic/Scharpery | Weight reduction critical |

### By Material Type

**Steel Chain:**
```julia
ph.materials["chain"] = PH.MaterialParameters(
    type="LinearElastic",
    E=2.1e11,
    μ=8.08e10
)
```
- **Pros**: Proven, robust, abrasion-resistant
- **Cons**: Heavy, expensive in deep water
- **Use**: Lower segments, harsh environments

**Wire Rope:**
```julia
ph.materials["wire"] = PH.MaterialParameters(
    type="LinearElastic",
    E=1.6e11,      # Lower than solid steel
    μ=6.15e10
)
```
- **Pros**: High strength-to-weight, flexible
- **Cons**: Fatigue-sensitive, corrosion
- **Use**: Intermediate segments, taut systems

**Polyester Rope:**
```julia
ph.materials["polyester"] = PH.MaterialParameters(
    type="LinearElastic",    # or "Scharpery" for time-dependent
    E=5.0e9,
    μ=1.92e9
)
```
- **Pros**: Lightweight, low cost, good fatigue
- **Cons**: Creep, axial stiffness variation
- **Use**: Upper segments, deep water

**Nylon Rope:**
```julia
ph.materials["nylon"] = PH.MaterialParameters(
    type="Scharpery",        # Significant viscoelasticity
    D0=4.0e-10,
    N=2,
    Dn=[2.0e-10, 1.0e-10],
    λn=[1.0e-2, 1.0e-3],
    ...
)
```
- **Pros**: High energy absorption, compliance
- **Cons**: Large creep, temperature-sensitive
- **Use**: Shock absorption, compliant systems

**HMPE (High-Modulus Polyethylene):**
```julia
ph.materials["hmpe"] = PH.MaterialParameters(
    type="LinearElastic",
    E=1.0e11,      # Very high for synthetic
    μ=3.85e10
)
```
- **Pros**: Very high strength-to-weight, low creep
- **Cons**: Expensive, UV-sensitive
- **Use**: Taut moorings, ultra-deep water

---

## Best Practices

### Material Property Selection

1. **Use realistic values**: Material properties significantly affect results
   ```julia
   # Bad: Using default values
   ph.materials["steel"] = PH.MaterialParameters(type="LinearElastic", E=1.0, μ=0.5)
   
   # Good: Using measured/literature values
   ph.materials["steel"] = PH.MaterialParameters(type="LinearElastic", E=2.1e11, μ=8.08e10)
   ```

2. **Account for submerged conditions**: Material properties can vary underwater
   ```julia
   # Synthetic ropes: properties may change with water absorption
   # Use wet properties, not dry properties
   E_dry = 6.0e9
   E_wet = 5.0e9  # Reduced due to water absorption
   
   ph.materials["polyester"] = PH.MaterialParameters(E=E_wet, ...)
   ```

3. **Consider temperature effects**: If applicable
   ```julia
   # For deep water applications, consider temperature
   # E typically decreases with increasing temperature
   ```

### Shear Modulus Calculation

Always compute shear modulus ($\mu$) consistently:

```julia
# Given Young's modulus and Poisson's ratio
E = 2.1e11      # Young's modulus [Pa]
ν = 0.3         # Poisson's ratio [-]

# Compute shear modulus
μ = E / (2 * (1 + ν))

ph.materials["material"] = PH.MaterialParameters(E=E, μ=μ, ...)
```

**Common Poisson's ratios:**
- Steel: ν ≈ 0.29-0.30
- Polyester: ν ≈ 0.35-0.40
- Nylon: ν ≈ 0.40-0.45

### Material Naming

Use descriptive, consistent tags:

```julia
# Good: Clear, descriptive
ph.materials["chain_76mm_grade3"] = ...
ph.materials["polyester_200mm"] = ...
ph.materials["wire_rope_spiral"] = ...

# Less clear
ph.materials["mat1"] = ...
ph.materials["material"] = ...
```

### Validation

Always validate material behavior:

```julia
# Check stress-strain response
E = 2.1e11
ε = 0.001  # 0.1% strain
σ_expected = E * ε

# Verify in post-processing that computed stresses are reasonable
```

### Viscoelastic Materials

When using Scharpery model:

1. **Start simple**: Use N=1 (single relaxation time) first
   ```julia
   ph.materials["poly_simple"] = PH.MaterialParameters(
       type="Scharpery",
       D0=2.0e-10,
       N=1,
       Dn=[1.0e-10],
       λn=[1.0e-3],
       ...
   )
   ```

2. **Calibrate from experiments**: Use creep or stress relaxation test data
3. **Consider time scales**: Choose relaxation times appropriate for your simulation
   ```julia
   # For quasi-static analysis (hours): λn ~ 1e-4 to 1e-2 [s⁻¹]
   # For dynamic analysis (minutes): λn ~ 1e-2 to 1e-1 [s⁻¹]
   ```

---

## See Also

- [API: Materials](../API/Physics/Materials.md)
- [User Guide: Segments](segments.md)
- [User Guide: Lines](lines.md)
- [API: TangentialDiffCalculus](../API/Physics/TangentialDiffCalculus.md)
