module gnlPara2D

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
using Interpolations
using WaveSpec.Constants
using Printf



"""
Warmup and Test params
===========

"""
@with_kw struct Warmup_params

  initCSV::String = "models/catShape_xfl60_zfl20.csv" 
  resDir = "data/"

  # Material properties
  E = 64.2986e9 #N
  mDry = 52.8  #kg/m Dry weight per unit len
  mSub = 45.936  #kg/m Submerged weight per unit len
  L = 75 #m
  ρcDry = 7.8e3 #kg/m3 Density of steel

  # Parameter Domain
  nx = 100

  # Time Parameters
  t0 = 0.0
  simT = 0.04
  simΔt = 0.02
  outΔt = 0.2

  # Fairlead Excitation
  fairLead_η = 0.1
  fairLead_T = 4.0

  # Drag coeff
  Cdn = 0.5
  Cdt = 0.5

end




"""
main()
======

"""
function main(params)

  @unpack initCSV, resDir = params
  pltName = resDir*"/gnl_"


  # Material properties
  @unpack E, ρcDry, L = params
  ρw = 1025 #Kg/m3 density of water
  μₘ = 0.5*E
  ρcSub = ρcDry - ρw
  @show L

  
  # Time Parameters
  @unpack t0, simT, simΔt, outΔt = params


  # Parameter domain
  @unpack nx = params
  domain = (0, L)
  partition = (nx)
  model = CartesianDiscreteModel(domain, partition)


  # Labelling 
  labels_Ω = get_face_labeling(model)
  add_tag_from_tags!(labels_Ω,"anchor",[1]) 
  add_tag_from_tags!(labels_Ω,"fairLead",[2]) 
  writevtk(model, pltName*"model")


  # Triangulations
  Ω = Interior(model)
  Γ = Boundary(model)


  ## Reference config Catenary
  # ---------------------Start---------------------    
  interpX, interpZ = setInitXZ(initCSV)

  X(r) = VectorValue( interpX(r[1]/L), interpZ(r[1]/L) )

  writevtk(Ω, pltName*"referenceDomain",
    cellfields=["X"=>X])
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


  ## Define Trial Fnc Static
  # ---------------------Start---------------------
  # Dirichlet BC
  gAnch_S(x) = VectorValue(0.0, 0.0)  

  US = TrialFESpace(Ψu, [gAnch_S, gAnch_S])
  # ----------------------End----------------------


  ## Define Trial Fnc Dynamic
  # ---------------------Start---------------------
  @unpack fairLead_η, fairLead_T = params

  # Dirichlet BC
  gAnch(x, t::Real) = VectorValue(0.0, 0.0)
  gAnch(t::Real) = x -> gAnch(x, t)
  gFairLead(x, t::Real) = VectorValue(fairLead_η*sin(2*pi/fairLead_T*t), 0.0)
  gFairLead(t::Real) = x -> gFairLead(x, t)

  U = TransientTrialFESpace(Ψu, [gAnch, gFairLead])
  # ----------------------End----------------------


  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓ = Measure(Γ,degree)
  @show nΓ = get_normal_vector(Γ)  


  # Initial solution
  Xh = interpolate_everywhere(X, Ψu)
  FWeih = interpolate_everywhere(VectorValue(0.0, -ρcSub*g), Ψu)

  
  ## Geometric quantities
  # ---------------------Start---------------------  
  J = ∇(Xh)' #d x q

  G = J' ⋅ J #q x q
  GInv = inv(G)

  Q = J ⋅ GInv

  P = (J ⋅ J')/ (J ⊙ J)  # Only true for q=1 d=2

  ∇X_Dir(u) = ∇(u)' ⋅ Q'

  FΓ(u) = ∇X_Dir(u) + TensorValue(1.0,0.0,0.0,1.0)

  JNew(u) = ∇(Xh + u)'
  sΛ(u) = ((JNew(u) ⊙ JNew(u)) ./ (J ⊙ J)).^0.5

  # Length of line
  calcLen(JLoc) = sum(∫( sqrt∘(JLoc ⊙ JLoc) )dΩ)

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
  
  
  ## Spring bed
  # ---------------------Start---------------------
  bedK1 = ρcSub*g
  bedK2 = ρcSub*g*1000
  bedRamp = 1e3
  spng(u) = 0.5+0.5*(tanh∘( bedRamp*excursion(u) ))
  excursion(u) = VectorValue(0.0,-1.0) ⋅ (Xh+u)
  bedForce(u) = spng(u) * (bedK1 + bedK2*excursion(u)) 
  # ----------------------End----------------------


  ## Drag force
  # ---------------------Start---------------------
  @unpack Cdn = params
  velN(v, u) = (v ⋅ t1(u)) * t1(u) - v
  velN_mag(v, u) = (velN(v, u) ⋅ velN(v, u)).^0.5
  # drag(v, u) = 0.5 * ρw * Cdn * velN(v, u) * velN_mag(v, u)  

  function drag(v, u)

    l_t1s = FΓ(u) ⋅ T1s
    l_t1 = l_t1s / ((l_t1s ⋅ l_t1s).^0.5)
  
    vn = (v ⋅ l_t1) * l_t1 - v
    vnm = (vn ⋅ vn).^0.5
    return 0.5 * ρw * Cdn * vn * vnm
  
  end
  # ----------------------End----------------------  


  ## Weak form
  # ---------------------Start---------------------
  mass(t, u, ψu) =  
  ∫( ( (ψu ⋅ ∂tt(u)) * ρcDry ) * ( (J ⊙ J).^0.5 ) )dΩ 


  damp(t, u, ψu) =  
  ∫( 
    (
      ψu ⋅ drag(∂t(u), u)
    )*((J ⊙ J).^0.5) 
  )dΩ 


  res(u, ψu) =  
  ∫( 
    ( 
      ∇X_Dir(ψu) ⊙ stressK(u) + 
      - ( ψu ⋅ FWeih ) + 
      - ( ψu ⋅ VectorValue(0.0,1.0) * bedForce(u) )
    )*((J ⊙ J).^0.5) 
  )dΩ 


  resD(t, u, ψu) = mass(t,u,ψu) + res(u,ψu) + damp(t, u, ψu)

  op_S = FEOperator(res, US, Ψu)
  op_D = TransientFEOperator(resD, U, Ψu; order=2)
  # ----------------------End----------------------


  ## Static Solution
  # ---------------------Start---------------------
  fx_U(r) = +0.00001*sin(π*r[1])
  fy_U(r) = -0.00001*sin(π*r[1])
  Ua(r) = VectorValue(fx_U(r), fy_U(r))  
  U0 = interpolate_everywhere(Ua, US)  

  nls = NLSolver(show_trace=true, 
    method=:newton, linesearch=Static(), 
    iterations=100, ftol = 1e-8, xtol = 1e-8)

  # nls = NLSolver(show_trace=true, 
  #   method=:anderson,  iterations=10)
  
  (uh_S, cache) = solve!(U0, nls, op_S)

  println("Length Catenary solution = ", calcLen(J) )  
  println("Length Static Solution = ", calcLen(JNew(uh_S)) )  

  xNew = Xh + uh_S

  writevtk(Ω, pltName*"staticRes",
    cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh_S, 
      "ETang"=>ETang(uh_S), "sigma"=>stressσ(uh_S),
      "t1"=>t1(uh_S) ])
  # ----------------------End----------------------  


  ## Initial solution
  # ---------------------Start---------------------  
  U0 = interpolate_everywhere(uh_S, U(t0))

  U0t = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))
  U0tt = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))
  # ----------------------End----------------------


  ## Dynamic Solver
  # ---------------------Start---------------------
  nls = NLSolver(show_trace=true, 
    method=:newton, linesearch=Static(), 
    iterations=100, ftol = 1e-8, xtol = 1e-8)

  ode_solver = GeneralizedAlpha(nls, simΔt, 0.0)    

  solnht = solve(ode_solver, op_D, (U0,U0t,U0tt), t0, simT) 
  # ----------------------End----------------------


  ## Save initial solution
  # ---------------------Start---------------------
  createpvd(pltName*"tSol") do pvd
    uh = U0
    tval = @sprintf("%5.6f",t0)                
    println("Time : $tval")
    tprt = @sprintf("%d",floor(Int64,t0*1000))                        
  
    xNew = X + uh
  
    pvd[t0] = createvtk(Ω,    
      pltName*"tSol_$tprt"*".vtu",
      cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh, 
        "ETang"=>ETang(uh), "sigma"=>stressσ(uh),
        "t1"=>t1(uh_S),
        "spr"=>spng(uh) ])
  end
  # ----------------------End----------------------


  ## Execute
  # ---------------------Start---------------------
  @show outMod = floor(Int64,outΔt/simΔt);
  tick()
  createpvd(pltName*"tSol", append=true) do pvd    
    cnt=0
    for (uh, t) in solnht                       
      cnt = cnt+1    
      tval = @sprintf("%5.6f",t)                
      println("Count \t $cnt \t Time : $tval")
      tprt = @sprintf("%d",floor(Int64,t*1000000))
            

      if(cnt%outMod == 0)               

        println("Paraview output")
      
        xNew = X + uh

        pvd[t] = createvtk(Ω,    
          pltName*"tSol_$tprt"*".vtu",
          cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh, 
            "ETang"=>ETang(uh), "sigma"=>stressσ(uh),
            "t1"=>t1(uh_S),
            "spr"=>spng(uh) ])
     
      end

      tock()
      println("-x-x-x-")
      println()
      tick()
    end
  end  
  tock()  
  # ----------------------End----------------------

    
end




"""
Aux functions
=============

"""
# ---------------------Start---------------------

function setInitXZ(initCSV)
  
  initXZ = CSV.read(initCSV, DataFrame, header=false)
  initXZ = Matrix(initXZ)

  dx = initXZ[2:end,:].-initXZ[1:end-1,:]
  dx = [0 0; dx]

  ds = zeros(size(dx,1))
  r = zeros(size(dx,1))

  for i in axes(dx,1)
    ds[i] = norm(dx[i,:])
  end

  for i in 2:length(ds)
    r[i] = r[i-1] + ds[i]
  end

  r = r / r[end]

  interpX = linear_interpolation(r, initXZ[:,1])
  interpZ = linear_interpolation(r, initXZ[:,2])
  
  # X(r) = VectorValue(interpX(r), interpZ(r))

  # Xh = interpolate_everywhere(X, Ψu)

  # return Xh

  return interpX, interpZ

end



# ramp = (0.0, 0.5/fFreq[2])
# rampΔT = ramp[2] - ramp[1]
# function getRamp(t::Real)
#   if(rampΔT < 1e-4)
#     return 1.0
#   end
  
#   if(t < ramp[1])
#     return 0.0
#   elseif(t > ramp[2])
#     return 1.0
#   else
#     return 0.5*( 1.0 - cos(π*(t - ramp[1]) / rampΔT) )
#   end
# end

# ----------------------End----------------------

end