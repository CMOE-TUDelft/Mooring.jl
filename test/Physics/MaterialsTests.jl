module MaterialsTests

import Mooring.Materials as M
import Mooring.TangentialDiffCalculus as TDC
using Gridap
using Test

# Test material properties
linear_elastic = M.LinearElastic(E=1.0)
@test linear_elastic.μ == 0.5

# Construct reference map CellField
model = CartesianDiscreteModel((0,1),(1,))
Ω = Interior(model)
dΩ = Measure(Ω,0)

# 2D map
reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},1)
V = TestFESpace(model,reffe)
X(r) = VectorValue(2r[1],0.0)
u(r) = VectorValue(r[1],0.0)
Xₕ = interpolate(X,V)
uₕ = interpolate(u,V)

# Define operators (2D)
Xq = get_cell_points(dΩ) 

# Test stresses
S = M.S(linear_elastic,Xₕ,uₕ)
FΓ = TDC.FΓ(uₕ,Xₕ)
K = M.K∘(FΓ,S)
J = TDC.J(Xₕ)
j = TDC.j(FΓ,J)
Λ = TDC.Λ(j,J)
σ = M.σ∘(Λ,FΓ,K)
@test S(Xq) == [[TensorValue{2,2,Float64}(0.625,0.0,0.0,0.0)]]
@test K(Xq) == [[TensorValue{2,2,Float64}(0.9375,0.0,0.0,0.0)]]
@test σ(Xq) == [[TensorValue{2,2,Float64}(0.9375,0.0,0.0,0.0)]]
end