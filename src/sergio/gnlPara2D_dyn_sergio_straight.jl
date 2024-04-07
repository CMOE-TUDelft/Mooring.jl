# module gnlPara2D

using DrWatson
@quickactivate "Mooring"


using Revise
using Gridap
using Gridap.Algebra
using DataFrames:DataFrame
using DataFrames:Matrix
using TickTock
using Parameters
using LineSearches: BackTracking
using LineSearches: Static
using LinearAlgebra
using CSV
using Printf
using Plots
using WaveSpec.WaveTimeSeries


function solveStatic()

  ## Define Trial Fnc
  # ---------------------Start---------------------
  # Dirichlet BC
  gSAnch(x) = gAnch(x, 0.0) 
  gSFairLead(x) = gFairLead(x, 0.0) 

  Us = TrialFESpace(Ψu, [gSAnch, gSFairLead])
  # ----------------------End----------------------

  op_Static = FEOperator(res, Us, Ψu)

  # Initial solution
  fx_U(r) = 0.01*sin(π*r[1])
  fy_U(r) = 0.0*sin(2*π*r[1])
  Ua(r) = VectorValue(fx_U(r), fy_U(r))
  U0 = interpolate_everywhere(Ua, Us)

  
  nls = NLSolver(show_trace=true, 
    method=:newton, linesearch=Static(), 
    iterations=100, ftol = 1e-8, xtol = 1e-8)


  (uh, cache) = solve!(U0, nls, op_Static)

  return uh
end


g = 9.81

# File properties
filename = datadir("sims_sergio","line_A0p10_f05p39_5","gnl")


# Material properties
E = 1e7  #N/m2
L = 10 #m
ρc = 100 #Kg/m3
μₘ = 0.5*E
@show g
@show E, ρc, L


# Line position
# xa = (0.0, 0.0)
# xf = (10.0, 0.0)


# Oscillation
fAmp = 0.1 #m
fFreq = 5.395 #Hz
fT = 1/fFreq #s
@show fAmp, fT


# Time Parameters
t0 = 0.0
simΔt = fT/20.0
simT = 80*fT
outΔt = fT/4.0
startRamp = (0.0, 4*fT)
@show startRamp


# Parameter Domain
nx = 100
domain = (0, L)
partition = (nx)
model = CartesianDiscreteModel(domain, partition)


# Labelling 
labels_Ω = get_face_labeling(model)
add_tag_from_tags!(labels_Ω,"anchor",[1]) 
add_tag_from_tags!(labels_Ω,"fairLead",[2]) 
writevtk(model, filename*"_model")


# Triangulations
Ω = Interior(model)
Γ = Boundary(model)


## Reference config 
# ---------------------Start---------------------
# sag = 0.10 #37.14594 #17.21725
# fx_X(r) = xa[1] + r[1]*xf[1]
# fy_X(r) = xa[2] + r[1]*xf[2] - sag*sin(π*r[1])

# # x^2
# fx_X(r) = xa[1] + r[1]*xf[1]
# fy_X(r) = xa[2] + r[1]*r[1]*xf[2]

# Straight
fx_X(r) = 0.0
fy_X(r) = r[1]


# # Fries paper 5.2
# fx_X(r) = r[1]
# fy_X(r) = 1/2*(1-r[1]) - 1/7*sin(π*r[1])

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
  dirichlet_tags=["anchor", "fairLead"])
# ----------------------End----------------------


## Define Trial Fnc
# ---------------------Start---------------------
# Dirichlet BC
gAnch(x, t::Real) = VectorValue(0.0, 0.0)
gAnch(t::Real) = x -> gAnch(x, t)

function gFairLead(x,t)    
  tRamp = timeRamp(t, startRamp[1], startRamp[2])
  fω = 2*pi*fFreq
  return VectorValue( tRamp*fAmp*sin(fω*t), 1.0 )
end   
gFairLead(t::Real) = x -> gFairLead(x, t)

U = TransientTrialFESpace(Ψu, [gAnch, gFairLead])
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

# Tangent
T1s = J ⋅ VectorValue(1.0)
T1m = (T1s ⋅ T1s).^0.5
T1 = T1s / T1m
t1s(u) = FΓ(u) ⋅ T1s
t1m(u) = (t1s(u) ⋅ t1s(u)).^0.5
t1(u) = t1s(u) / t1m(u)

# Strain tensors
EDir(u) = 0.5 * ( FΓ(u)' ⋅ FΓ(u) - TensorValue(1.0,0.0,0.0,1.0) )
ETang(u) = P ⋅ EDir(u) ⋅ P

# Stress tensors
stressK(u) = 2*μₘ * (FΓ(u) ⋅ ETang(u))
stressS(u) = 2*μₘ * ETang(u)
stressσ(u) = ( FΓ(u) ⋅ stressS(u) ⋅ FΓ(u)' ) / sΛ(u)
# ----------------------End----------------------


## Cell state
# ---------------------Start---------------------  
function new_J(J_csin, Jin)
  return true, Jin
end

function create_cellState(Jin, loc)
  local J_cs
  J_cs = CellState(Jin(loc), dΩ)  
  update_state!(new_J, J_cs, Jin)    

  return J_cs
end

loc = Point(0.0)
Xh_cs = create_cellState(Xh, loc)
J_cs = create_cellState(J, loc)
QTrans_cs = create_cellState(Q', loc)
P_cs = create_cellState(P, loc)
JJ_cs = create_cellState((J ⊙ J).^0.5, loc)  
T1s_cs = create_cellState(T1s, loc)
T1m_cs = create_cellState(T1m, loc)
T1_cs = create_cellState(T1, loc)
# ----------------------End----------------------


## Function form stressK
# ---------------------Start---------------------
# stressK function for speed
function stressK_fnc(u)
  
  local FΓ, EDir, ETang
  
  FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)

  EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )

  ETang = P_cs ⋅ EDir ⋅ P_cs

  return 2*μₘ * (FΓ ⋅ ETang)

end
# ----------------------End----------------------  


# # Body force
# FBodyh = interpolate_everywhere(VectorValue(0.0, 0.0), Ψu)


## Weak form
# ---------------------Start---------------------
mass(t, u, ψu) =  
  ∫( ( (ψu ⋅ ∂tt(u)) * ρc ) * JJ_cs )dΩ 

res(u, ψu) =  
  ∫( ( (∇(ψu)' ⋅ QTrans_cs) ⊙ stressK_fnc(u) )*JJ_cs )dΩ  


resD(t, u, ψu) = mass(t,u,ψu) + res(u,ψu)


op_AD = TransientFEOperator(resD, U, Ψu; order=2)
# ----------------------End----------------------


## Solver
# ---------------------Start---------------------

nls = NLSolver(show_trace=true, 
  method=:newton, linesearch=Static(), 
  iterations=100, ftol = 1e-8, xtol = 1e-8)

ode_solver = GeneralizedAlpha(nls, simΔt, 0.0)    
# ----------------------End----------------------


## Initial solution
# ---------------------Start---------------------
println()
println("-----Static Solution-----")
u0 = solveStatic()
# u0 = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))

LNew = sum(∫( sqrt∘(JNew(u0) ⊙ JNew(u0)) )dΩ)
println("Length Original = ", LOrig )  
println("Length New = ", LNew )  
println("Ratio = ", LNew/LOrig ) 
println("-----End Static Solution-----")
println()

# da = wload(datadir("sims","mor_sol.jld2"))
# uf = FEFunction( Ψu, da["uh"] )
# U0 = interpolate_everywhere(uf, U(t0))

u0t = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))
u0tt = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))
# ----------------------End----------------------


solnht = solve(ode_solver, op_AD, (u0,u0t,u0tt), t0, simT) 


## Probing
# ---------------------Start---------------------
nPrb = 11
xPrb = LinRange(domain[1], domain[2], nPrb)
xPrb = Point.(xPrb)

daFile1 = open( filename*"_data1.dat", "w" )
# ----------------------End----------------------


function getPrincipalStress(uh, lprb)
  sT = stressσ(uh).(lprb)
  sTP = [ (lT[1] + lT[4])/2 + sqrt( ((lT[1] - lT[4])/2)^2 + lT[2]^2 )
    for lT in sT ] 
  return sTP
end

function getPrincipalETang(uh, lprb)  
  sETang = ETang(uh).(lprb)  
  sETangP = [ (lT[1] + lT[4])/2 + sqrt( ((lT[1] - lT[4])/2)^2 + lT[2]^2 )
    for lT in sETang ] 
  return sETangP
end


# Save initial soln
createpvd(filename*"_tSol") do pvd
  uh = u0
  tval = @sprintf("%5.6f",t0)                
  println("Time : $tval")
  tprt = @sprintf("%d",floor(Int64,t0*1000))                        

  xNew = X + uh

  xNew_prb = xNew.(xPrb)  
  sT_prb = getPrincipalStress(uh, xPrb)

  pvd[t0] = createvtk(Ω,    
      filename*"_tSol_$tprt"*".vtu",
      cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh, 
      "ETang"=>ETang(uh), "sigma"=>stressσ(uh)])
end




# Execute
tick()
@show outMod = floor(Int64,outΔt/simΔt);
createpvd(filename*"_tSol", append=true) do pvd    
  cnt=0
  for (uh, t) in solnht                       
    cnt = cnt+1        
    waveCycle = t/fT
    @printf("Time : %10.4f s \t Counter : %5i \t Wave : %5.2f \n", t, cnt, waveCycle)    
    tprt = @sprintf("%d",floor(Int64,t*1000000))                        

    xNew = X + uh

    xNew_prb = xNew.(xPrb)  
    sETang_prb = getPrincipalETang(uh, xPrb)
    sT_prb = getPrincipalStress(uh, xPrb)

    @printf(daFile1, "%15.6f", t)
    
    [@printf(daFile1, ", %15.6f, %15.6f, %20.10e, %20.10e",
        lx[1], lx[2], lETang, lsT) 
      for (lx, lETang, lsT) in zip(xNew_prb, sETang_prb, sT_prb)]

    @printf(daFile1, "\n")
    
    println("-x-x-")
    tock()
    println()    
    tick()    

    if(cnt%outMod != 0) 
      continue
    end    

    pvd[t] = createvtk(Ω,    
      filename*"_tSol_$tprt"*".vtu",
      cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh, 
      "ETang"=>ETang(uh), "sigma"=>stressσ(uh)])
      
  end
end  
tock()  
close(daFile1)


daF = CSV.read(filename*"_data1.dat", DataFrame, header=false)
daF = Matrix(daF)



function setGridlines(plt)
  vline!(plt, [startRamp[2]/fT], linewidth = 3, lc = "red", label=nothing)
  plot!(plt, grid=:true, gridcolor=:black, 
    gridalpha=0.5, gridlinestyle=:dot,
    gridlinewidth = 1)    
  plot!( plt, dpi = 330,
    titlefontsize=15,
    tickfontsize=13, 
    labelfontsize=15,
    legendfontsize = 13)    
  plot!( plt, 
    xlabel = "t/T")
end


prbi = round(Int64, (nPrb+1)/2)
plt1 = plot()
plot!( plt1, daF[:,1]/fT, daF[:,1+(prbi-1)*4+1], label = nothing, linewidth=3 )
plot!( plt1, title = "X Coord", ylabel = "(m)" )
setGridlines(plt1)

prbi = round(Int64, (nPrb+1)/2)
plt2 = plot()
plot!( plt2, daF[:,1]/fT, daF[:,1+(prbi-1)*4+2], label = nothing, linewidth=3 )
plot!( plt2, title = "Z Coord", ylabel = "(m)" )
setGridlines(plt2)

prbi = 1
plt3 = plot()
plot!( plt3, daF[:,1]/fT, daF[:,1+(prbi-1)*4+3], label = nothing, linewidth=3 )
plot!( plt3, title = "ETang at Anchor", ylabel = "ETang" )
setGridlines(plt3)

prbi = nPrb
plt4 = plot()
plot!( plt4, daF[:,1]/fT, daF[:,1+(prbi-1)*4+3], label = nothing, linewidth=3 )
plot!( plt4, title = "ETang at FairLead", ylabel = "ETang" )
setGridlines(plt4)


prbi1 = 1
prbi2 = nPrb
plt5 = plot()
plot!( plt5, daF[:,1]/fT, daF[:,1+(prbi2-1)*4+3]-daF[:,1+(prbi1-1)*4+3], label = nothing, linewidth=3 )
plot!( plt5, title = "ETang FairLead-Anchor", ylabel = "ETang" )
setGridlines(plt5)

savefig(plt1, filename*"_plot_posX.png")
savefig(plt2, filename*"_plot_posZ.png")
savefig(plt3, filename*"_plot_ETang_Anchor.png")
savefig(plt4, filename*"_plot_ETang_FairLead.png")
savefig(plt5, filename*"_plot_ETang_FairLeadVsAnchor.png")


# end