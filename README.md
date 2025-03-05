# Mooring.jl

This repository is built using the [Julia Language](https://julialang.org/) and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to facilitate a reproducible scientific project named:

> **Mooring.jl**

Authored by **Shagun**.

## Getting Started

To reproduce this project locally, follow these steps:

### 1. Clone the Repository
First, download the codebase by cloning this repository:

```sh
git clone https://github.com/yourusername/Mooring.jl.git
cd Mooring.jl
```

**Note:** Raw data files are typically **not** included in version control (git). You may need to download them separately.

### 2. Start a Julia Session
Open a terminal in the project directory and start Julia with:

```sh
julia --project=.
```

### 3. Install Dependencies
Once in Julia, enter package manager mode by pressing `]` and run:

```julia
(Mooring) pkg> instantiate
```

This will install all required dependencies, ensuring that everything runs correctly, including proper resolution of local paths.

---

## Updating Specific Dependencies
If you need to update specific dependencies to newer versions, enter package manager mode (`]`) and run:

```julia
(Mooring) pkg> add https://github.com/shagunTUD/Gridap.jl#devMoor
(Mooring) pkg> add https://github.com/shagun751/WaveSpec.jl#main
```

This will fetch and install the specified versions from their respective repositories.

---

## Running Scripts
You can execute a script using:

```julia
julia> include("test/demo_test.jl")
```

### DrWatson Integration
Most scripts in this project begin with:

```julia
using DrWatson
@quickactivate "Mooring.jl"
```

This automatically activates the project environment and ensures proper handling of local paths.

---

## Additional Notes
- If you run into issues with missing dependencies, try running `] resolve` in Julia to fix package compatibility.
- Ensure that your Julia version matches the one specified in `Project.toml` to avoid compatibility issues.
