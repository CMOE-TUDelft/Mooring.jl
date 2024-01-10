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
using Printf


g = 9.81

# File properties
filename = datadir("sims","run","mem")


# Material properties
E = 1e4
ρc = 2000/g #Kg/m3
L = 1 #m
A = 0.01
μₘ = 0.5*E
@show EA = E*A
@show m = ρc * A
@show L

# Time Parameters
t0 = 0.0
simΔt = 0.02
simT = 5.0
outΔt = 0.02


# Parameter Domain
nx = 100
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
g1(x, t::Real) = VectorValue(0.0, 0.0)
g1(t::Real) = x -> g1(x, t)
g2(x, t::Real) = VectorValue(-0.02*sin(2*pi/0.2*t), 0.0)
g2(t::Real) = x -> g2(x, t)

U = TransientTrialFESpace(Ψu, [g2,g1])
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


FBodyh = interpolate_everywhere(VectorValue(0.0, -ρc*g), Ψu)

bedK = ρc*g*2
bedRamp = 1e3
spng(u) = 0.5+0.5*(tanh∘( VectorValue(0.0,-bedRamp) ⋅ (Xh+u)))

## Weak form
# ---------------------Start---------------------
mass(t, u, ψu) =  
  ∫( ( (ψu ⋅ ∂tt(u)) * ρc ) * ( (J ⊙ J).^0.5 ) )dΩ 

res(u, ψu) =  
  ∫( 
    ( 
      ∇X_Dir(ψu) ⊙ stressK(u) + 
      - ( ψu ⋅ FBodyh ) + 
      - ( ψu ⋅ VectorValue(0.0,bedK) * spng(u) )
    )*((J ⊙ J).^0.5) 
  )dΩ 


resD(t, u, ψu) = mass(t,u,ψu) + res(u,ψu)


op_AD = TransientFEOperator(resD, U, Ψu; order=2)
# ----------------------End----------------------


## Initial solution
# ---------------------Start---------------------

# Initial solution
da = wload(datadir("sims","mem_sol.jld2"))
uf = FEFunction( Ψu, da["uh"] )

# Ua(r) = VectorValue(0.0, -0.001*sin(π*r[1]))
# U0 = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))

U0 = interpolate_everywhere(uf, U(t0))

U0t = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))
U0tt = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))
# ----------------------End----------------------


## Solver
# ---------------------Start---------------------

nls = NLSolver(show_trace=true, 
  method=:newton, linesearch=Static(), 
  iterations=100, ftol = 1e-8)

ode_solver = GeneralizedAlpha(nls, simΔt, 0.0)    

solnht = solve(ode_solver, op_AD, (U0,U0t,U0tt), t0, simT) 
# ----------------------End----------------------


createpvd(filename*"_tSol") do pvd
  uh = U0
  tval = @sprintf("%5.6f",t0)                
  println("Time : $tval")
  tprt = @sprintf("%d",floor(Int64,t0*1000))                        

  xNew = X + uh

  pvd[t0] = createvtk(Ω,    
      filename*"_tSol_$tprt"*".vtu",
      cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh, 
      "ETang"=>ETang(uh), "sigma"=>stressσ(uh),
      "spr"=>spng(uh) ])
end



# Execute
tick()
@show outMod = floor(Int64,outΔt/simΔt);
createpvd(filename*"_tSol", append=true) do pvd    
  cnt=0
  for (uh, t) in solnht                       
    cnt = cnt+1    
    tval = @sprintf("%5.6f",t)                
    println("Time : $tval")
    tprt = @sprintf("%d",floor(Int64,t*1000000))                        

    if(cnt%outMod != 0) 
      continue
    end
    xNew = X + uh

    pvd[t] = createvtk(Ω,    
      filename*"_tSol_$tprt"*".vtu",
      cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh, 
      "ETang"=>ETang(uh), "sigma"=>stressσ(uh),
      "spr"=>spng(uh) ])
  end
end  
tock()  

# end