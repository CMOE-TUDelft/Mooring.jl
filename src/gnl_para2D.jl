# module gnlPara2D

using DrWatson
@quickactivate "Mooring"


using Revise
using Gridap
using Gridap.Algebra
using Plots
using DataFrames:DataFrame
using DataFrames:Matrix
using TickTock
using Parameters
using LineSearches: BackTracking
using LineSearches: Static
using LinearAlgebra
using CSV


# File properties
filename = datadir("sims","mem")


# Material properties
E = 1e4
ρₘ = 1
L = 1 #m
A = 0.01
μₘ = 0.5*E
@show EA = E*A
@show m = ρₘ * A
@show L


# Parameter Domain
nx = 1000
domain = (0, L)
partition = (nx)
model = CartesianDiscreteModel(domain, partition)


# Labelling 
labels_Ω = get_face_labeling(model)
add_tag_from_tags!(labels_Ω,"leftEdge",[1]) 
add_tag_from_tags!(labels_Ω,"rightEdge",[2]) 
writevtk(model, filename*"_model")


# Triangulations
Ω = Interior(model)
Γ = Boundary(model)


# ## Reference config
# # ---------------------Start---------------------
# @show X_sag = 0.01*L
# @show αX1 = 0.6
# αX2 = (X_sag*(2/αX1/L)^2, αX1*L/2, -X_sag)
# @show αX2
# fx_X(r) = αX1*r[1]
# fy_X(r) = αX2[1] * (fx_X(r)-αX2[2]).^2 + αX2[3]
# X(r) = VectorValue(fx_X(r), fy_X(r))

# writevtk(Ω, filename*"_referenceDomain",
#   cellfields=["X"=>X])
# # ----------------------End----------------------


## Reference config Fries paper 5.2
# ---------------------Start---------------------
fx_X(r) = r[1]
fy_X(r) = 1/2*(1-r[1]) - 1/7*sin(π*r[1])
# fy_X(r) = -r[1]
X(r) = VectorValue(fx_X(r), fy_X(r))

# writevtk(Ω, filename*"_referenceDomain",
#   cellfields=["X"=>X])
# ----------------------End----------------------


## Define Test Fnc
# ---------------------Start---------------------
order = 1

reffe = ReferenceFE(lagrangian, 
  VectorValue{2,Float64}, order)
Ψu = TestFESpace(Ω, reffe, 
  conformity=:H1, 
  dirichlet_tags=["leftEdge", "rightEdge"])
# ----------------------End----------------------


## Define Trial Fnc
# ---------------------Start---------------------
# Dirichlet BC
g1(x) = VectorValue(0.0, 0.0)

U = TrialFESpace(Ψu, [g1,g1])
# ----------------------End----------------------


# Measures
degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
@show nΓ = get_normal_vector(Γ)


## Geometric quantities
# ---------------------Start---------------------
Xh = interpolate_everywhere(X, Ψu)

J = ∇(Xh)'

LOrig = sum(∫( sqrt∘(J ⊙ J) )dΩ)

println("Reference Length = ", LOrig )  

G = J' ⋅ J
GInv = inv(G)

Q = J ⋅ GInv

P = (J ⋅ J')/ (J ⊙ J)  # Only true for q=1 d=2

∇X_Dir(u) = ∇(u)' ⋅ Q'

FΓ(u) = ∇X_Dir(u) + TensorValue(1.0,0.0,0.0,1.0)

JNew(u) = ∇(Xh + u)'
sΛ(u) = ((JNew(u) ⊙ JNew(u)) ./ (J ⊙ J)).^0.5

# Strain tensors
EDir(u) = 0.5 * ( FΓ(u)' ⋅ FΓ(u) - TensorValue(1.0,0.0,0.0,1.0) )
ETang(u) = P ⋅ EDir(u) ⋅ P

# Stress tensors
stressK(u) = 2*μₘ * (FΓ(u) ⋅ ETang(u))
stressS(u) = 2*μₘ * ETang(u)
stressσ(u) = ( FΓ(u) ⋅ stressS(u) ⋅ FΓ(u)' ) / sΛ(u)

# ----------------------End----------------------


FExth = interpolate_everywhere(VectorValue(0.0, -2000), Ψu)


## Weak form
# ---------------------Start---------------------
res(u, ψu) =  
  ∫( ( ∇X_Dir(ψu) ⊙ stressK(u) - ( ψu ⋅ FExth ) )*((J ⊙ J).^0.5) )dΩ 

# ----------------------End----------------------


## Solver
# ---------------------Start---------------------

# op = FEOperator(res, U, Ψu)
# uh = solve(op)

# # Initial solution
# @show U_sag = 0.01*L
# αU2 = (U_sag*(2/αX1/L)^2, αX1*L/2, -U_sag)
# @show αU2
# fx_U(r) = 0.0
# fy_U(r) = αU2[1] * (fx_X(r)-αX2[2]).^2 + αU2[3]
# Ua(r) = VectorValue(fx_U(r), fy_U(r))
# # Ua(r) = VectorValue(0.0, -0.00001)
# U0 = interpolate_everywhere(Ua, U)

# Initial solution
fx_U(r) = 0.0#sin(π*r[1])
fy_U(r) = -0.01*sin(π*r[1])
Ua(r) = VectorValue(fx_U(r), fy_U(r))
# Ua(r) = VectorValue(0.0, -0.00001)
U0 = interpolate_everywhere(Ua, U)

op = FEOperator(res, U, Ψu)

# nls = NewtonRaphsonSolver(LUSolver(), 1e-5, 100)

nls = NLSolver(show_trace=true, 
  method=:newton, linesearch=Static(), 
  iterations=100, ftol = 1e-8)

# nls = NLSolver(show_trace=true, 
#   method=:anderson,  iterations=10)

# solver = FESolver(nls)

#uh = solve(nls, op)
(uh, cache) = solve!(U0, nls, op)


xNew = Xh + uh

writevtk(Ω, filename*"_referenceDomain",
  cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh, 
    "ETang"=>ETang(uh), "sigma"=>stressσ(uh) ])


LNew = sum(∫( sqrt∘(JNew(uh) ⊙ JNew(uh)) )dΩ)
println("Length Original = ", LOrig )  
println("Length New = ", LNew )  
println("Ratio = ", LNew/LOrig )  

# Potential energy calculation
pot = A/2 * sum(∫( ETang(uh) ⊙ stressS(uh) * ((J ⊙ J).^0.5) )dΩ)
println("Stored potential energy = ", pot)

# ----------------------End----------------------

# end