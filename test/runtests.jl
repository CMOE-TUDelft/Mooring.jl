using Test

# Run test suite
println("Starting tests")

@time @testset "EnvironmentalConditions" begin include(joinpath("EnvironmentalConditionsTests.jl")) end
@time @testset "SeaBed" begin include(joinpath("SeaBedTests.jl")) end
@time @testset "Drag" begin include(joinpath("DragTests.jl")) end
@time @testset "PointMotion" begin include(joinpath("PointMotionTests.jl")) end
# @time @testset "Demo" begin include(joinpath("demo_test.jl")) end