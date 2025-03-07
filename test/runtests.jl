using Test

# Run test suite
println("Starting tests")

@time @testset "EnvironmentalConditions" begin include(joinpath("EnvironmentalConditionsTests.jl")) end
@time @testset "SeaBed" begin include(joinpath("SeaBedTests.jl")) end
# @time @testset "Demo" begin include(joinpath("demo_test.jl")) end