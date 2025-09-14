# Installation

To install the latest version of Mooring.jl:

## Using `Pkg` REPL mode (recommended)

1. Open Julia  
2. Enter the package manager by typing `]`  
3. Add the package with:

```julia
pkg> add https://github.com/CMOE-TUDelft/Mooring.jl
```

## Using `Pkg` in scripts

Alternatively, install directly from the Julia REPL or in your scripts:

```julia
using Pkg
Pkg.add(url="https://github.com/CMOE-TUDelft/Mooring.jl")
```

## Developing locally

If you want to work on the code and documentation:

```bash
git clone https://github.com/CMOE-TUDelft/Mooring.jl
cd Mooring.jl
julia --project=docs/ -e 'using Pkg; Pkg.instantiate()'
```

