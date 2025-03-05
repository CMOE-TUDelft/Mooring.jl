using DrWatson, Test
@quickactivate "Mooring"

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
println("Starting tests")

@time @testset "Demo" begin include(projectdir("test","demo_test.jl")) end