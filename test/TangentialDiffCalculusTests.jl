module TangentialDiffCalculusTests

import Mooring.TangentialDiffCalculus as TDC
using Gridap
using Test

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

# Test operators (2D)
Xq = get_cell_points(dΩ) # Quadrature point at the center of the element (0.5,)
J = TDC.J(Xₕ)
G = TDC.G(J)
Q = TDC.Q(J)
∇ₓΓdir = TDC.∇ₓΓdir(uₕ,Xₕ)
@test J(Xq) == [[TensorValue{2,1,Float64}(2.0,0.0)]]
@test G(Xq) == [[TensorValue{1,1,Float64}(4.0)]]
@test Q(Xq) == [[TensorValue{2,1,Float64}(0.5,0.0)]]
@test ∇ₓΓdir(Xq) == [[TensorValue{2,2,Float64}(0.5,0.0,0.0,0.0)]]

# 3D map
reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},1)
V = TestFESpace(model,reffe)
X(r) = VectorValue(2r[1],0.0,0.0)
u(r) = VectorValue(r[1],0.0,0.0)
Xₕ = interpolate(X,V)
uₕ = interpolate(u,V)

# Test operators (3D)
J = TDC.J(Xₕ)
G = TDC.G(J)
Q = TDC.Q(J)
∇ₓΓdir = TDC.∇ₓΓdir(uₕ,Xₕ)
@test J(Xq) == [[TensorValue{3,1,Float64}(2.0,0.0,0.0)]]
@test G(Xq) == [[TensorValue{1,1,Float64}(4.0)]]
@test Q(Xq) == [[TensorValue{3,1,Float64}(0.5,0.0,0.0)]]
@test ∇ₓΓdir(Xq) == [[TensorValue{3,3,Float64}(0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)]]

end