# Mooring.jl
<img src="docs/src/assets/logo.svg" width="250" title="Mooring.jl logo">

**Mooring.jl** is a Julia package for simulating mooring lines using the Finite Element Method. The package is based on **dynamic finite strain theory** and **tangential differential calculus**. **Mooring.jl** supports linear elastic, nonlinear elastic and nonlinear viscoelastic material properties, providing an accurate and efficient framework for modeling the behavior of mooring lines in offshore applications, either for chain/steel based cables or synthetic ropes.

## Features
- **Finite Element Method (FEM)** for high-fidelity mooring line simulations, with arbitrary order of interpolation
- **Dynamic finite strain theory** for large deformation analysis
- **Tangential differential calculus** for geometric consistency in simulations with complex material models
- **Linear and nonlinear material models** including viscoelastic effects, relevant for synthetic ropes
- **Efficient time integration schemes** for dynamic analysis
- **Modular and extensible** design for research and engineering applications

## Installation

You can install **Mooring.jl** using Julia's package manager:

```julia
using Pkg
Pkg.add("Mooring")
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

Contributions are welcome! If you find a bug or have a feature request, please open an issue or submit a pull request. Please, use the [CONTRIBUTING.md](CONTRIBUTING.md) guidelines.

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
