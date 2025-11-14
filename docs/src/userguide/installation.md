# Installation

This guide covers installing Mooring.jl and its dependencies.

## Requirements

**Mooring.jl** requires:
- **Julia 1.9 or higher** (recommended: Julia 1.10+)
- Operating systems: Linux, macOS, or Windows
- Sufficient RAM: minimum 4GB, recommended 8GB+ for large problems

## Standard Installation

### Using `Pkg` REPL mode (recommended)

This is the simplest method for most users:

1. **Open Julia**  
   Launch Julia from your terminal or application launcher

2. **Enter package manager**  
   Type `]` to enter the Pkg REPL mode (the prompt changes to `pkg>`)

3. **Add the package**
   ```julia
   pkg> add https://github.com/CMOE-TUDelft/Mooring.jl
   ```

4. **Exit package mode**  
   Press backspace or `Ctrl+C` to return to the Julia REPL

5. **Verify installation**
   ```julia
   julia> using Mooring
   julia> println("Mooring.jl installed successfully!")
   ```

### Using `Pkg` in scripts

If you prefer scripting or need programmatic installation:

```julia
using Pkg
Pkg.add(url="https://github.com/CMOE-TUDelft/Mooring.jl")
```

This method is useful for:
- Automated setup scripts
- Containerized environments
- CI/CD pipelines

## Dependencies

Mooring.jl automatically installs its dependencies, including:

### Core Dependencies
- **[Gridap.jl](https://github.com/gridap/Gridap.jl)**: Finite element framework
- **[GridapGmsh.jl](https://github.com/gridap/GridapGmsh.jl)**: Gmsh mesh support
- **[Parameters.jl](https://github.com/mauro3/Parameters.jl)**: Parameter management
- **[Roots.jl](https://github.com/JuliaMath/Roots.jl)**: Root finding for catenary calculations

### Special Dependencies
Some features require custom branches:
- **WaveSpec.jl**: For wave spectrum analysis
- **Modified Gridap.jl**: With mooring-specific extensions

These are automatically handled during installation.

## Development Installation

For contributors or users who want to modify the source code:

### Clone the repository

```bash
git clone https://github.com/CMOE-TUDelft/Mooring.jl
cd Mooring.jl
```

### Set up the development environment

```julia
# In Julia REPL
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### Build documentation locally

```bash
julia --project=docs/ -e 'using Pkg; Pkg.instantiate()'
julia --project=docs/ docs/make.jl
```

The built documentation will be in `docs/build/`.

### Run tests

```julia
# In Julia REPL from project root
using Pkg
Pkg.test("Mooring")
```

Or from the command line:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Troubleshooting

### Package resolution issues

If you encounter dependency conflicts:

```julia
using Pkg
Pkg.resolve()
Pkg.instantiate()
```

### Precompilation errors

If precompilation fails, try:

```julia
using Pkg
Pkg.build("Mooring")
Pkg.precompile()
```

### Updating the package

To update to the latest version:

```julia
using Pkg
Pkg.update("Mooring")
```

### Clean reinstall

If problems persist, perform a clean reinstall:

```julia
using Pkg
Pkg.rm("Mooring")
Pkg.gc()  # Clean up unused packages
Pkg.add(url="https://github.com/CMOE-TUDelft/Mooring.jl")
```

## Verifying Installation

Test your installation with a simple example:

```julia
using Mooring
import Mooring.ParameterHandlers as PH

# Create a parameter handler
ph = PH.ParameterHandler()

# If this runs without errors, installation was successful!
println("✓ Mooring.jl is ready to use")
```

## Next Steps

After successful installation:
1. Review the [Quick Start Example](../UserGuide.md#quick-start-example)
2. Learn about [Points](points.md) in the mooring system
3. Explore the [API Documentation](../API/API.md)

## Getting Help

If you encounter installation issues:
- Check the [GitHub Issues](https://github.com/CMOE-TUDelft/Mooring.jl/issues)
- Ask questions in [GitHub Discussions](https://github.com/CMOE-TUDelft/Mooring.jl/discussions)
- Review Julia's [package management documentation](https://docs.julialang.org/en/v1/stdlib/Pkg/)

