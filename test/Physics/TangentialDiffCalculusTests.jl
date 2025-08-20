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

# Define operators (2D)
Xq = get_cell_points(dΩ) # Quadrature point at the center of the element (0.5,)
J = TDC.J(Xₕ)
G = TDC.G(J)
G_composed = TDC.G∘(J)
Q = TDC.Q(J)
Q_composed = TDC.Q∘(J)
∇ₓΓdir = TDC.∇ₓΓdir(uₕ,Xₕ)
∇ₓΓdir_composed = TDC.∇ₓΓdir∘(∇(uₕ),J)
FΓ = TDC.FΓ(uₕ,Xₕ)
FΓ_composed = TDC.FΓ∘(∇ₓΓdir)
P = TDC.P(J)
P_composed = TDC.P∘(J)
j = TDC.j(FΓ,J)
Λ = TDC.Λ(j,J)
Λ_composed = TDC.Λ∘(j,J)
Edir = TDC.Edir∘(FΓ)
Etang = TDC.Etang∘(P,Edir)

# Test operators (2D)
@test J(Xq) == [[TensorValue{2,1,Float64}(2.0,0.0)]]
@test G(Xq) == [[TensorValue{1,1,Float64}(4.0)]]
@test G_composed(Xq) == [[TensorValue{1,1,Float64}(4.0)]]
@test Q(Xq) == [[TensorValue{2,1,Float64}(0.5,0.0)]]
@test Q_composed(Xq) == [[TensorValue{2,1,Float64}(0.5,0.0)]]
@test ∇ₓΓdir(Xq) == [[TensorValue{2,2,Float64}(0.5,0.0,0.0,0.0)]]
@test ∇ₓΓdir_composed(Xq) == [[TensorValue{2,2,Float64}(0.5,0.0,0.0,0.0)]]
@test FΓ(Xq) == [[TensorValue{2,2,Float64}(1.5,0.0,0.0,1.0)]]
@test FΓ_composed(Xq) == [[TensorValue{2,2,Float64}(1.5,0.0,0.0,1.0)]]
@test P(Xq) == [[TensorValue{2,2,Float64}(1.0,0.0,0.0,0.0)]]
@test P_composed(Xq) == [[TensorValue{2,2,Float64}(1.0,0.0,0.0,0.0)]]
@test Λ(Xq) == [[1.5]]
@test Λ_composed(Xq) == [[1.5]]
@test Edir(Xq) == [[TensorValue{2,2,Float64}(0.625,0.0,0.0,0.0)]]
@test Etang(Xq) == [[TensorValue{2,2,Float64}(0.625,0.0,0.0,0.0)]]

# 3D map
reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},1)
V = TestFESpace(model,reffe)
X(r) = VectorValue(2r[1],0.0,0.0)
u(r) = VectorValue(r[1],0.0,0.0)
Xₕ = interpolate(X,V)
uₕ = interpolate(u,V)

# Define operators (3D)
J = TDC.J(Xₕ)
G = TDC.G(J)
G_composed = TDC.G∘(J)
Q = TDC.Q∘(J)
Q_composed = TDC.Q(J)
∇ₓΓdir = TDC.∇ₓΓdir(uₕ,Xₕ)
∇ₓΓdir_composed = TDC.∇ₓΓdir∘(∇(uₕ),J)
FΓ = TDC.FΓ(uₕ,Xₕ)
FΓ_composed = TDC.FΓ∘(∇ₓΓdir)
P = TDC.P(J)
P_composed = TDC.P∘(J)
j = TDC.j(FΓ,J)
Λ = TDC.Λ(j,J)
Λ_composed = TDC.Λ∘(j,J)
Edir = TDC.Edir∘(FΓ)
Etang = TDC.Etang∘(P,Edir)

# Test operators (3D)
@test J(Xq) == [[TensorValue{3,1,Float64}(2.0,0.0,0.0)]]
@test G(Xq) == [[TensorValue{1,1,Float64}(4.0)]]
@test G_composed(Xq) == [[TensorValue{1,1,Float64}(4.0)]]
@test Q(Xq) == [[TensorValue{3,1,Float64}(0.5,0.0,0.0)]]
@test Q_composed(Xq) == [[TensorValue{3,1,Float64}(0.5,0.0,0.0)]]
@test ∇ₓΓdir(Xq) == [[TensorValue{3,3,Float64}(0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)]]
@test ∇ₓΓdir_composed(Xq) == [[TensorValue{3,3,Float64}(0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)]]
@test FΓ(Xq) == [[TensorValue{3,3,Float64}(1.5,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0)]]
@test FΓ_composed(Xq) == [[TensorValue{3,3,Float64}(1.5,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0)]]
@test P(Xq) == [[TensorValue{3,3,Float64}(1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)]]
@test P_composed(Xq) == [[TensorValue{3,3,Float64}(1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)]]
@test Λ(Xq) == [[1.5]]
@test Λ_composed(Xq) == [[1.5]]
@test Edir(Xq) == [[TensorValue{3,3,Float64}(0.625,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)]]
@test Etang(Xq) == [[TensorValue{3,3,Float64}(0.625,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)]]

end