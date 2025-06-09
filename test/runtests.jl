using Test

# Run test suite
println("Starting tests")

# Physics tests
@time @testset "EnvironmentalConditions" begin include(joinpath("Physics/EnvironmentalConditionsTests.jl")) end
@time @testset "SeaBed" begin include(joinpath("Physics/SeaBedTests.jl")) end
@time @testset "Drag" begin include(joinpath("Physics/DragTests.jl")) end
@time @testset "PointMotion" begin include(joinpath("Physics/PointMotionTests.jl")) end
@time @testset "TangentialDiffCalculus" begin include(joinpath("Physics/TangentialDiffCalculusTests.jl")) end
@time @testset "Materials" begin include(joinpath("Physics/MaterialsTests.jl")) end
# @time @testset "Demo" begin include(joinpath("demo_test.jl")) end