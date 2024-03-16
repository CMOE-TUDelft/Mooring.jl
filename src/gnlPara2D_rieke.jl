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
using Printf
using WaveSpec.Constants
using WaveSpec.Jonswap
using WaveSpec.WaveTimeSeries
using WaveSpec.Currents



"""
Warmup and Test params
===========

"""
@with_kw struct Test_params

  initCSV::String = "models/catShape_xfl60_zfl20.csv" 
  resDir = "data/"

  # Material properties
  E = 64.2986e9 #N
  L = 75 #m
  A_str = 2*(π*0.048*0.048/4) #m2 Str cross-section area
  ρcDry = 7.8e3 #kg/m3 Density of steel  

  # Bed parameter
  bed_tanhRamp = 1e3
  bed_springK = 0.5 #30kN/m4
  #based on Marco's suggestion of 30kN/m2/m in the OrcaFlex manual

  # Parameter Domain
  nx = 100
  order = 1

  # Time Parameters
  t0 = 0.0
  simT = 0.04
  simΔt = 0.02
  outΔt = 0.2
  maxIter = 100

  # Drag coeff
  C_dn = 2.6 # Normal drag coeff
  d_dn = 0.048 #m Normal drag projection diameter
  C_dt = 1.4 # Tangent drag coff
  d_dt = 0.048/pi #m Tangent drag projection diameter

  # Added mass coeff
  C_an = 1.0 # Normal added-mass coeff
  C_at = 0.5 # Tangent added-mass coff

  # Time signal ramp up (t0 t1)
  startRamp = (0.0, 15)

  # Wave spectrum
  h0 = 23 #m
  Hs = 3 #m
  Tp = 12 #s
  nω = 257 #including 0

  # Current
  strCur = CurrentStat(23, [-23.0, -11.0, 0.0], [0.0, 0.0, 0.0])

end




"""
main()
======

"""
function main(params)

  @unpack initCSV, resDir = params
  pltName = resDir*"/gnl_"


  # Material properties
  @unpack E, ρcDry, L, A_str = params
  ρw = 1025 #Kg/m3 density of water
  μₘ = 0.5*E
  ρcSub = ρcDry - ρw
  @show L

  
  # Time Parameters
  @unpack t0, simT, simΔt, outΔt, maxIter = params


  ## Wave input
  # ---------------------Start---------------------  
  @unpack Hs, Tp, h0, nω = params
  ω, S, A = jonswap(Hs, Tp,
    plotflag=false, nω = nω)

  k = dispersionRelAng.(h0, ω; msg=false)
  α = randomPhase(ω)

  sp = SpecStruct( h0, ω, S, A, k, α; Hs = Hs, Tp = Tp )
  #η, ϕ, u, w = waveAiry1D(sp, t, 0.1, -0.1)
  # ----------------------End----------------------  


  ## Mesh setup: Regular
  # ---------------------Start---------------------  
  @unpack nx = params
  domain = (0, L)
  partition = (nx)      
  model = CartesianDiscreteModel(domain, partition)

  # Labelling 
  labels_Ω = get_face_labeling(model)
  add_tag_from_tags!(labels_Ω,"anchor",[1]) 
  add_tag_from_tags!(labels_Ω,"fairLead",[2]) 
  writevtk(model, pltName*"model")

  daFile0 = open( pltName*"aout.dat", "w" )
  # ----------------------End----------------------  


  # ## Mesh setup: Gmesh
  # # ---------------------Start---------------------  
  # model = DiscreteModelFromFile(projectdir("models/mesh2/irrMesh.msh"))

  # # Labelling 
  # labels_Ω = get_face_labeling(model)
  # # add_tag_from_tags!(labels_Ω,"anchor",[1]) 
  # # add_tag_from_tags!(labels_Ω,"fairLead",[2]) 
  # writevtk(model, pltName*"model")
  # # ----------------------End----------------------  


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
  @unpack order = params

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
  gFairLead_S(x) = VectorValue(0.0, 0.0)  

  US = TrialFESpace(Ψu, [gAnch_S, gFairLead_S])
  # ----------------------End----------------------


  ## Define Trial Fnc Dynamic
  # ---------------------Start---------------------  
  # Dirichlet BC
  gAnch(x, t::Real) = VectorValue(0.0, 0.0)
  gAnch(t::Real) = x -> gAnch(x, t)

  @show Xh_fl = X(Point(L))
  @show (Xh_fl[1], Xh_fl[2]-h0)
  @unpack startRamp = params
  function getFairLeadEnd(x,t)    
    η, px, py = waveAiry1D_pPos(sp, t, Xh_fl[1], Xh_fl[2]-h0)
    tRamp = timeRamp(t, startRamp[1], startRamp[2])

    # return VectorValue(0.0, η*tRamp)
    return VectorValue(px*tRamp, py*tRamp)
  end
  
  # function getFairLeadEnd(x,t)    
  #   fairLead_η = 1.5 #m2
  #   fairLead_ω = 1 #rad/s

  #   tRamp = timeRamp(t, startRamp[1], startRamp[2])

  #   return VectorValue( -tRamp*fairLead_η*sin(fairLead_ω*t), 
  #   tRamp*fairLead_η*(1-cos(fairLead_ω*t)) )
  # end    
  
  gFairLead(x, t::Real) = getFairLeadEnd(x,t)    
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
  
  
  function getPrincipalStress(uh, lp)

    sT = stressσ(uh).(lp)
    sTP = [ (lT[1] + lT[4])/2 + sqrt( ((lT[1] - lT[4])/2)^2 + lT[2]^2 )
      for lT in sT ]     
  
  end
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
  FWeih_cs = create_cellState(FWeih, loc)
  J_cs = create_cellState(J, loc)
  QTrans_cs = create_cellState(Q', loc)
  P_cs = create_cellState(P, loc)
  JJ_cs = create_cellState((J ⊙ J).^0.5, loc)  
  T1s_cs = create_cellState(T1s, loc)
  T1m_cs = create_cellState(T1m, loc)
  T1_cs = create_cellState(T1, loc)
  # ----------------------End----------------------
  
  
  ## Spring bed
  # ---------------------Start---------------------
  @unpack bed_tanhRamp, bed_springK = params
  @show bed_tanhRamp, bed_springK
  bedK1 = ρcSub*g
  bedK2 = ρcSub*g * bed_springK
  println("Bed spring constant = ", bedK2, " N/m3/m")
  bedRamp = bed_tanhRamp
  spng(u) = 0.5+0.5*(tanh∘( bedRamp*excursion(u) ))
  excursion(u) = VectorValue(0.0,-1.0) ⋅ (Xh+u)
  bedForce(u) = spng(u) * (bedK1 + bedK2*excursion(u)) 

  function bedSpring_fnc(u)
    local exc, lspng

    exc = VectorValue(0.0,-1.0) ⋅ (Xh_cs+u)
    lspng = 0.5 + 0.5*(tanh∘( bedRamp*exc ))

    return lspng * (bedK1 + bedK2*exc) 
  end

  function bedDamp_fnc(v, u)
    local exc, lspng, vz

    exc = VectorValue(0.0,-1.0) ⋅ (Xh_cs + u)
    lspng = 0.5 + 0.5*(tanh∘( bedRamp*exc ))

    vz = VectorValue(0.0,1.0) ⋅ (v)

    return -lspng * (0.05*bedK2*vz)
  end
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


  ## Function form drag
  # ---------------------Start---------------------
  @unpack C_dn, d_dn, C_dt, d_dt = params
  @unpack strCur = params
  @show D_dn = 0.5 * ρw * C_dn * d_dn / A_str #kg/m4
  @show D_dt = 0.5 * ρw * C_dt * π * d_dt / A_str #kg/m4

  # Constant Currents
  function getCurrentField(r)
    pz = Xh(r) ⋅ VectorValue(0.0,1.0) - strCur.h0
    pz = min(pz, 0.0)
    pz = max(pz, -strCur.h0)

    return VectorValue( strCur.itp( pz ), 0.0 ) 
  end
  UCur_h = interpolate_everywhere(getCurrentField, Ψu)
  UCur_cs = create_cellState( UCur_h, loc )

  function drag_ΓX_intp(v,u)  #No wave and current
    # Slow basic version
    # Useful for plotting

    local vn, vnm, vt, vtm
    
    vt = -(v ⋅ t1(u)) * t1(u)
    vtm = (vt ⋅ vt).^0.5
    vn = -v - vt
    vnm = (vn ⋅ vn).^0.5    

    return (D_dn * vn * vnm + D_dt * vt * vtm) * sΛ(u) 

  end

  function drag_ΓX(v, u) #No wave and current

    local FΓ, t1s, t1m2, vn, vnm, sΛ, vt, vtm

    FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)

    t1s = FΓ ⋅ T1s_cs
    t1m2 = t1s ⋅ t1s    
    #t1 = t1s / ((t1s ⋅ t1s).^0.5)

    sΛ = (t1m2.^0.5) / T1m_cs
    
    vt = -(v ⋅ t1s) * t1s / t1m2
    vtm = (vt ⋅ vt).^0.5
    vn = -v - vt
    vnm = (vn ⋅ vn).^0.5    

    return (D_dn * vn * vnm + D_dt * vt * vtm) * sΛ 

  end

  function drag_n_ΓX(t, v, u)

    local FΓ, t1s, t1m2, vn, vnm, sΛ, vr, tRamp

    tRamp = timeRamp(t, startRamp[1], startRamp[2])

    FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)

    t1s = FΓ ⋅ T1s_cs
    t1m2 = t1s ⋅ t1s    
    #t1 = t1s / ((t1s ⋅ t1s).^0.5)

    sΛ = (t1m2.^0.5) / T1m_cs
    
    vr = UCur_cs * tRamp - v
    vn = vr - (vr ⋅ t1s) * t1s / t1m2
    # vn = (v ⋅ t1s) * t1s / t1m2 - v 
    vnm = (vn ⋅ vn).^0.5

    return D_dn * vn * vnm * sΛ 

  end

  function drag_t_ΓX(t, v, u)

    local FΓ, t1s, t1m2, sΛ, vt, vtm, tRamp

    tRamp = timeRamp(t, startRamp[1], startRamp[2])

    FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)

    t1s = FΓ ⋅ T1s_cs
    t1m2 = t1s ⋅ t1s    
    #t1 = t1s / ((t1s ⋅ t1s).^0.5)

    sΛ = (t1m2.^0.5) / T1m_cs
    
    vt = ((UCur_cs * tRamp - v) ⋅ t1s) * t1s / t1m2
    # vt = -(v ⋅ t1s) * t1s / t1m2    
    vtm = (vt ⋅ vt).^0.5

    return D_dt * vt * vtm * sΛ 

  end
  # ----------------------End----------------------  


  ## Function form added-mass
  # ---------------------Start---------------------
  @unpack C_an, C_at, = params
  @show D_an = ρw * C_an
  @show D_at = ρw * C_at   

  function addedMass_n_ΓX(a, u)

    local FΓ, t1s, t1m2, an, sΛ

    FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)

    t1s = FΓ ⋅ T1s_cs
    t1m2 = t1s ⋅ t1s    
    #t1 = t1s / ((t1s ⋅ t1s).^0.5)

    sΛ = (t1m2.^0.5) / T1m_cs
    
    an = (a ⋅ t1s) * t1s / t1m2 - a    

    return D_an * an * sΛ

  end

  function addedMass_t_ΓX(a, u)

    local FΓ, t1s, t1m2, sΛ, at

    FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)

    t1s = FΓ ⋅ T1s_cs
    t1m2 = t1s ⋅ t1s    
    #t1 = t1s / ((t1s ⋅ t1s).^0.5)

    sΛ = (t1m2.^0.5) / T1m_cs
    
    at = -(a ⋅ t1s) * t1s / t1m2
    
    return D_at * at * sΛ 

  end
  # ----------------------End----------------------  


  ## Weak form: Static
  # ---------------------Start---------------------
  # res(u, ψu) =  
  # ∫( 
  #   ( 
  #     ∇X_Dir(ψu) ⊙ stressK_fnc(u) + 
  #     - ( ψu ⋅ FWeih_cs ) + 
  #     - ( ψu ⋅ VectorValue(0.0,1.0) * bedForce(u) )
  #   )*JJ_cs 
  # )dΩ 
  
  res(u, ψu) =      
    ∫( ( (∇(ψu)' ⋅ QTrans_cs) ⊙ stressK_fnc(u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ FWeih_cs )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ VectorValue(0.0,1.0) * bedSpring_fnc(u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ drag_n_ΓX(0.0, VectorValue(0.0,0.0), u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ drag_t_ΓX(0.0, VectorValue(0.0,0.0), u) )*JJ_cs )dΩ 


  op_S = FEOperator(res, US, Ψu)
  # ----------------------End----------------------


  ## Weak form: Dynamic
  # ---------------------Start---------------------
  # resD(t, u, ψu) =  
  # ∫( 
  #   ( 
  #     (ψu ⋅ ∂tt(u)) * ρcDry +
  #     ∇X_Dir(ψu) ⊙ stressK(u) +
  #     -ψu ⋅ FWeih +
  #     -ψu ⋅ VectorValue(0.0,1.0) * bedForce(u) 
  #   )*((J ⊙ J).^0.5)
  # )dΩ 


  resD(t, u, ψu) =  
    ∫( ( (ψu ⋅ ∂tt(u)) * ρcDry )*JJ_cs )dΩ +
    # ∫( ( -ψu ⋅ addedMass_n_ΓX(∂tt(u), u) )*JJ_cs )dΩ +
    # ∫( ( -ψu ⋅ addedMass_t_ΓX(∂tt(u), u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ VectorValue(0.0,1.0) * bedDamp_fnc(∂t(u), u) )*JJ_cs )dΩ +
    ∫( ( (∇(ψu)' ⋅ QTrans_cs) ⊙ stressK_fnc(u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ FWeih_cs )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ VectorValue(0.0,1.0) * bedSpring_fnc(u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ drag_n_ΓX(t, ∂t(u), u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ drag_t_ΓX(t, ∂t(u), u) )*JJ_cs )dΩ 
    # ∫( ( -ψu ⋅ drag_ΓX(∂t(u), u) )*JJ_cs )dΩ 
    # ∫( (  )*JJ_cs )dΩ +
    # ∫( (  )*JJ_cs )dΩ +
    # ∫( (  )*JJ_cs )dΩ   


  op_D = TransientFEOperator(resD, U, Ψu; order=2)
  # ----------------------End----------------------


  ## Static Solution
  # ---------------------Start---------------------
  fx_U(r) = +0.00001*sin(π*r[1])
  fy_U(r) = -0.00001*sin(π*r[1])
  Ua(r) = VectorValue(fx_U(r), fy_U(r))  
  U0 = interpolate_everywhere(Ua, US)    
  
  # nls = NLSolver(LUSolver(), show_trace=true, 
  nls = NLSolver(show_trace=true, 
    method=:newton, linesearch=Static(), 
    iterations=maxIter, ftol = 1e-8, xtol = 1e-8)

  # nls = NLSolver(show_trace=true, 
  #   method=:anderson,  iterations=10)
  
  (uh_S, cache) = solve!(U0, nls, op_S)

  println("Length Catenary solution = ", calcLen(J) )  
  println("Length Static Solution = ", calcLen(JNew(uh_S)) )  

  xNew = Xh + uh_S

  writevtk(Ω, pltName*"staticRes",
    cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh_S, 
      "ETang"=>ETang(uh_S), "sigma"=>stressσ(uh_S) ])
  # ----------------------End----------------------  


  ## Initial solution
  # ---------------------Start---------------------  
  U0 = interpolate_everywhere(uh_S, U(t0))

  U0t = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))
  U0tt = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))
  # ----------------------End----------------------


  ## Dynamic Solver
  # ---------------------Start---------------------
  # nls = NLSolver(LUSolver(), show_trace=true, 
  nls = NLSolver(show_trace=true, 
    method=:newton, linesearch=Static(), 
    iterations=maxIter, ftol = 1e-8, xtol = 1e-8)

  # nls = NewtonRaphsonSolver(LUSolver(), 1e-8, 100)

  ode_solver = GeneralizedAlpha(nls, simΔt, 0.0)    

  solnht = solve(ode_solver, op_D, (U0,U0t,U0tt), t0, simT) 
  # ----------------------End----------------------


  ## Save quantities
  # ---------------------Start---------------------  
  rPrb = Point.(0.0:L/10:L)
  nPrb = length(rPrb)

  daFile1 = open( pltName*"data1.dat", "w" )
  # daFile2 = open( pltName*"data2.dat", "w" )
  
  daSave1 = DataFrame(zeros(Float64, 1, 9), :auto)
  daSave2 = DataFrame(zeros(Float64, 1, 6), :auto)
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
        "spr"=>spng(uh) ])
        # "drag_ΓX" => drag_ΓX_intp(uh,uh) ])
  end
  # ----------------------End----------------------


  ## Execute
  # ---------------------Start---------------------
  @show outMod = floor(Int64,outΔt/simΔt);
  execTime = zeros(Float64, 1, 10)
  execTime[1] = time()  # tick()
  execTime[3] = time()   
  tick()
  createpvd(pltName*"tSol", append=true) do pvd    
    cnt=0
    for (uh, t) in solnht                       
      cnt = cnt+1    
      tval = @sprintf("%5.6f",t)                
      println("Count \t $cnt \t Time : $tval")
      tprt = @sprintf("%d",floor(Int64,t*1000000))

      xNew = X + uh
      σT = getPrincipalStress(uh, rPrb)      
      xNewPrb = xNew.(rPrb)

      # lDa = [t, xPrb[1][1]
      #   xNew(xPrb[1])[1], xNew(xPrb[1])[2], σT[1], σT[1]*A_str,
      #   xNew(xPrb[end])[1], xNew(xPrb[end])[2], σT[end], σT[end]*A_str]        
      # push!(daSave1, lDa)

      @printf(daFile1, "%.5f \t",t)
      # [print(daFile1, string(val)*", \t") for val in lDa]
      [@printf(daFile1, "%.5f \t %.5f \t %.5f \t %.5f \t", 
        rPrb[i][1], xNewPrb[i][1], xNewPrb[i][2], σT[i])
        for i in 1:nPrb]
      @printf(daFile1, "\n")

            
      if(cnt%outMod == 0)               

        println("Paraview output")        

        pvd[t] = createvtk(Ω,    
          pltName*"tSol_$tprt"*".vtu",
          cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh, 
            "ETang"=>ETang(uh), "sigma"=>stressσ(uh),
            "spr"=>spng(uh) ])
            # "drag_ΓX" => drag_ΓX_intp(uh,uh) ])
     
      end
      execTime[4] = time()  
      tock()
      @printf(daFile0, 
        "Step Time: \t %5i \t %.3f \n", 
        cnt, execTime[4]-execTime[3])
      println("-x-x-x-")
      println()
      execTime[3] = time()  
      tick()
    end
  end  
  execTime[2] = time()  
  tock()
  @printf(daFile0, 
    "\nTotal Time: \t %5i \t %.3f \n", 
    round(simT/simΔt), execTime[2]-execTime[1])

  close(daFile0)
  close(daFile1)
  # close(daFile2)
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
# ----------------------End----------------------

end