# Mooring.jl

**Mooring.jl** is a Julia package for simulating mooring lines using the Finite Element Method. The package is based on **dynamic finite strain theory** and **tangential differential calculus**. **Mooring.jl** supports linear elastic, nonlinear elastic and nonlinear viscoelastic material properties, providing an accurate and efficient framework for modeling the behavior of mooring lines in offshore applications, either for chain/steel based cables or synthetic ropes.

#### Table of Contents

- [Key Features](#key-features)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)
- [Authors](#authors)
- [Contact](#contact)
- [Known Issues](#known-issues)

## Key Features

### 🔬 Advanced Finite Element Method
- **Arbitrary interpolation order**: Use linear, quadratic, or higher-order finite elements for enhanced accuracy
- **Wide variety of FE spaces**: Support for H¹-conforming spaces (continuous), L²-conforming spaces (discontinuous), and other Sobolev spaces
- **Flexible discretization**: Independent mesh refinement and interpolation order per segment
- **Multifield formulation**: Each mooring line uses multifield FE spaces combining segment-level spaces
- **Independent FE space per segment**: Modular design enabling different discretizations for each segment
- **Gridap.jl foundation**: Leverages powerful FEM abstractions for complex PDEs

### ⚙️ Robust Nonlinear Mechanics
- **Dynamic finite strain theory**: Geometrically exact formulation for arbitrarily large deformations
- **Nonlinear constitutive models**: Support for linear elastic and nonlinear viscoelastic materials (Scharpery model). Easy to extend to user-defined models
- **Tangential differential calculus**: Rigorous geometric framework ensuring consistency in 1D cable dynamics in 2D/3D space
- **Fully implicit time integration**: Robust schemes for highly nonlinear and dynamic problems with strong stability. Wide variety of schemes: Newmark-β, generalized-α, backward Euler, Runge-Kutta, ...
- **Automatic differentiation**: Jacobian computation for both geometric and material nonlinearities without manual derivations

### 🧩 Modular and Extensible Architecture
- **Easy material implementation**: Add custom constitutive models by implementing stress functions
- **Flexible forcing terms**: Simple to incorporate add-hoc forcing functions. Currently, drag forces, seabed interaction, buoyancy, wave loading are implemented
- **Clear separation**: Parameter configuration (user level) vs. FEM operations (solver level)

### 💻 User-Friendly Interface
- **YAML configuration**: Define complete mooring systems without writing code
- **Parameter handler**: Intuitive API for programmatic setup in Julia
- **Comprehensive documentation**: Detailed user guides and API reference
- **Easy installation**: Simple package manager installation
- **Minimal dependencies**: Core Julia ecosystem with Gridap.jl for FEM

## Installation

To install the latest version of **Mooring.jl**:

### Using `Pkg` REPL mode (recommended)

1. Open Julia  
2. Enter the package manager by typing `]`  
3. Add the package with:

```julia
pkg> add https://github.com/CMOE-TUDelft/Mooring.jl
```

### Using `Pkg` in scripts

Alternatively, install directly from the Julia REPL or in your scripts:

```julia
using Pkg
Pkg.add(url="https://github.com/CMOE-TUDelft/Mooring.jl")
```

## Usage

### Basic Example

```julia
using Mooring

# TO DO
```

## Results

### Example Simulation Output

Below is an example of a mooring line under wave loading:

*(Include plots or numerical results here)*

## Contributing

Contributions are welcome! If you find a bug or have a feature request, please open an issue or submit a pull request. Please, use the [CONTRIBUTING.md](https://github.com/CMOE-TUDelft/Mooring.jl/CONTRIBUTING.md) guidelines.

## License

This package is licensed under the MIT License.

## Authors

This repository has been developed by [Shagun Agarwal](https://github.com/shagun751) and [Oriol Colomés](https://github.com/oriolcg). See also the list of [Contributors](https://github.com/CMOE-TUDelft/Mooring.jl/graphs/contributors).

## Contact

For questions or collaboration, please reach out to the authors or post an [issue](https://github.com/CMOE-TUDelft/Mooring.jl/issues) to the repository.

---

## Known issues

### Updating Specific Dependencies
If you need to update specific dependencies to newer versions, enter package manager mode (`]`) and run:

```julia
(Mooring) pkg> add https://github.com/shagunTUD/Gridap.jl#devMoor
(Mooring) pkg> add https://github.com/shagun751/WaveSpec.jl#main
```

This will fetch and install the specified versions from their respective repositories.

### Additional Notes
- If you run into issues with missing dependencies, try running `] resolve` in Julia to fix package compatibility.
- Ensure that your Julia version matches the one specified in `Project.toml` to avoid compatibility issues.
