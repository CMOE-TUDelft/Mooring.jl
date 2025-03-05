using Test

# Run test suite
println("Starting tests")

@time @testset "Demo" begin include(joinpath("demo_test.jl")) end