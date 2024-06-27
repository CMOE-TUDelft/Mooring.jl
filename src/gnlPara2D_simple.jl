module gnlPara2D

using DrWatson
@quickactivate "Mooring"


using Revise
using Gridap
using Gridap.Algebra
using Gridap.ODEs
using Gridap.Arrays: testitem, return_cache
using DataFrames:DataFrame
using DataFrames:Matrix
using WriteVTK
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
  bed_springK = 50 #30kN/m4
  # Based on Marco's suggestion of 30kN/m2/m in the OrcaFlex manual
  # Following certain tests from Rieke database,
  # I choose to keep this value as bed_springK = 50 
  # => actual spring constant of 3.32 MN/m3/m

  # Parameter Domain
  nx = 100
  order = 1

  outFreeSurface = false

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
  seed = -1
  ωc = -1

  # Current
  strCur = CurrentStat(23, [-23.0, -11.0, 0.0], [0.0, 0.0, 0.0])

  # Forced fairlead motion
  ffm_η = 0.1 #m
  ffm_ω = 0.5 #Hz

end




"""
main()
======

"""
function main(params)

  @unpack initCSV, resDir = params
  pltName = resDir*"/gnl_"


  ## Terminal Output
  # ---------------------Start---------------------  
  outFile0 = open( pltName*"aout.dat", "w" )
  outFile1 = open( pltName*"arunTime.dat", "w" )
  
  printTer(a::String,b) = printTerAndFile(a,b,outFile0)
  printTer(a::String) = printTerAndFile(a,outFile0)
  printTer() = printTerAndFile("",outFile0)
  # ----------------------End----------------------  


  # Material properties
  @unpack E, ρcDry, L, A_str = params
  ρw = 1025 #Kg/m3 density of water
  μₘ = 0.5*E
  ρcSub = ρcDry - ρw
  printTer("[VAL] Length = ", L)
  printTer()
  
  # Time Parameters
  @unpack t0, simT, simΔt, outΔt, maxIter, startRamp = params
  @unpack outFreeSurface = params
  outMod = floor(Int64,outΔt/simΔt);
  
  printTer("[VAL] (t0, simT) = ",(t0, simT))
  printTer("[VAL] (simΔt, outΔt) = ",(simΔt, outΔt))
  printTer("[VAL] StartRamp = ",startRamp)
  printTer("[VAL] outMod = ", outMod)
  printTer()

  ## Wave input
  # ---------------------Start---------------------  
  @unpack h0 = params

  sp = getInputSpec(params)
  #η, ϕ, u, w = waveAiry1D(sp, t, 0.1, -0.1)

  printTer("[VAL] h0 = ", h0)
  printTer("[VAL] Hs = ", sp.Hs)
  printTer("[VAL] Tp = ", sp.Tp)
  printTer("[VAL] nw = ", sp.nω)
  printTer("[VAL] wc = ", sp.ω[end])
  printTer()
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

  printTer("[VAL] nx = ", nx)
  printTer()
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


  ## Domain to output surface elevation in Paraview
  # ---------------------Start---------------------
  model_fs = CartesianDiscreteModel((0, 1.2*L), (nx))
  Ω_fs = Interior(model_fs)

  Ψu_fs = FESpace(Ω_fs, 
    ReferenceFE(lagrangian, Float64, 1), 
    conformity=:H1)  
  
  X_fs(r) = r[1]
  
  function getEta_fs(t,x)    
    lx = x ⋅ VectorValue(1.0)
    return waveAiry1D_eta(sp, t, lx, 0.0)         
  end
  getEta_fs(t) = x -> getEta_fs(t,x)
  # ----------------------End----------------------  


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

  Xh_fl = X(Point(L))
  printTer("[VAL] Xh_fl = ", (Xh_fl[1], Xh_fl[2]))
  printTer("[VAL] Xh_fl = ", (Xh_fl[1], Xh_fl[2]-h0))
  printTer()

  # Create interpolable obj for this
  function createInterpObj(sp, t, x, z)
    η, px, py = waveAiry1D_pPos(sp, t, x, z)    

    tRamp = timeRamp.(t, startRamp[1], startRamp[2])

    itp = Interpolations.interpolate(
      η.*tRamp, 
      BSpline(Cubic(Line(OnGrid()))) )
    sitp_η = scale(itp, t)

    itp = Interpolations.interpolate(
      px.*tRamp, 
      BSpline(Cubic(Line(OnGrid()))) )
    sitp_px = scale(itp, t)

    itp = Interpolations.interpolate(
      py.*tRamp, 
      BSpline(Cubic(Line(OnGrid()))) )
    sitp_py = scale(itp, t)

    return sitp_η, sitp_px, sitp_py
  end

  sitp_η, sitp_px, sitp_py = 
    createInterpObj( sp, 
      t0:simΔt/2.0:(1.2*simT) , 
      Xh_fl[1], Xh_fl[2]-h0 )

  # function getFairLeadEnd(x,t)    

  #   # tRamp already done in createInterpObj
  #   px = sitp_px(t)
  #   py = sitp_py(t)

  #   # return VectorValue(0.0, η)
  #   return VectorValue(px, py)
  # end

  # function getFairLeadEnd(x,t)    
  #   η, px, py = waveAiry1D_pPos(sp, t, Xh_fl[1], Xh_fl[2]-h0)
  #   tRamp = timeRamp(t, startRamp[1], startRamp[2])

  #   # return VectorValue(0.0, η*tRamp)
  #   return VectorValue(px*tRamp, py*tRamp)
  # end
  
  @unpack ffm_η, ffm_ω = params
  printTer("[VAL] ffm_η", ffm_η)
  printTer("[VAL] ffm_ω", ffm_ω)
  function getFairLeadEnd(x,t)    
    # ffm_η = 1.5 #m2
    # ffm_ω = 1 #rad/s

    tRamp = timeRamp(t, startRamp[1], startRamp[2])    

    return VectorValue( 
      tRamp*ffm_η*sin(ffm_ω*t),
      0.0 )
  end    
  
  gFairLead(x, t::Real) = getFairLeadEnd(x,t)    
  gFairLead(t::Real) = x -> gFairLead(x, t)

  U = TransientTrialFESpace(Ψu, [gAnch, gFairLead])
  # ----------------------End----------------------


  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓ = Measure(Γ,degree)
  nΓ = get_normal_vector(Γ)  


  # Initial solution
  Xh = interpolate_everywhere(X, Ψu)
  FWeih = interpolate_everywhere(VectorValue(0.0, -ρcSub*g), Ψu)
  
  ## Geometric quantities
  # ---------------------Start---------------------  
  J = ∇(Xh)' #d x q

  G = J' ⋅ J #q x q
  GInv = inv(G)

  Q = J ⋅ GInv
  QTrans = Q'

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

  function getPrincipalStress2(lT)
    
    return (lT[1] + lT[4])/2 + sqrt( ((lT[1] - lT[4])/2)^2 + lT[2]^2 )
  
  end
  # ----------------------End----------------------


  ## Cell state
  # ---------------------Start---------------------  
  # function new_J(J_csin, Jin)
  #   return true, Jin
  # end

  function create_cellState(Jin, loc)
    local J_cs
    J_cs = CellState(Jin(loc), dΩ)  
    update_state!( (a,b) -> (true, b), J_cs, Jin)    

    return J_cs
  end

  loc = Point(0.0)
  Xh_cs = create_cellState(Xh, loc)
  Xh_x_cs = create_cellState(Xh⋅VectorValue(1.0,0.0), loc)
  Xh_z_cs = create_cellState(Xh⋅VectorValue(0.0,1.0), loc)
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
  bedK1 = ρcSub*g
  bedK2 = ρcSub*g * bed_springK
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

  printTer("[VAL] bed_tanhRamp", bed_tanhRamp)
  printTer("[VAL] bed_springK", bed_springK)
  printTer("[VAL] Bed spring constant (N/m3/m) = ", bedK2)
  printTer()
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

  function stressK_fnc(QTr, P, ∇u)
    
    local FΓ, EDir, ETang
    
    FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)

    # FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)

    EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )

    ETang = P ⋅ EDir ⋅ P

    return 2*μₘ * (FΓ ⋅ ETang)

  end


  function stressσ_fnc(QTr, P, J, ∇u)      
    
    local FΓ, EDir, ETang, stressS, JNew, sΛ
    
    FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)

    # FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)

    EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )

    ETang = P ⋅ EDir ⋅ P

    stressS = 2*μₘ * ETang

    JNew = J  + ∇u'
    sΛ = ((JNew ⊙ JNew) ./ (J ⊙ J)).^0.5

    return ( FΓ ⋅ stressS ⋅ FΓ' ) / sΛ

  end
  # ----------------------End----------------------  


  ## Wave velocity vector CellField
  # ---------------------Start---------------------
  function getWaveVel(r, t, xNew_loc)

    xqp = xNew_loc(r)
    tRamp = timeRamp(t, startRamp[1], startRamp[2])

    w_u, w_w = waveAiry1D_vel(sp, t, 
      xqp[1], xqp[2]-sp.h0 )

    return VectorValue(w_u, w_w) * tRamp
  end

  getWaveVel_cf(t, xNew_loc) = 
    CellField( r -> getWaveVel(r, t, xNew_loc), Ω )  
  waveVel_cs = 
    create_cellState( getWaveVel_cf(t0, Xh), loc)  
  # ----------------------End----------------------  


  ## Function form drag
  # ---------------------Start---------------------
  @unpack C_dn, d_dn, C_dt, d_dt = params
  @unpack strCur = params
  D_dn = 0.5 * ρw * C_dn * d_dn / A_str #kg/m4
  D_dt = 0.5 * ρw * C_dt * π * d_dt / A_str #kg/m4
  
  printTer("[VAL] D_dn = ", D_dn)
  printTer("[VAL] D_dt = ", D_dt)
  printTer()

  # Constant Currents
  function getCurrentField(r)
    pz = Xh(r) ⋅ VectorValue(0.0,1.0) - strCur.h0
    pz = min(pz, 0.0)
    pz = max(pz, -strCur.h0)

    return VectorValue( strCur.itp( pz ), 0.0 ) 
  end
  UCur_cs = create_cellState( 
    CellField( r -> getCurrentField(r), Ω), 
    loc )

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
    # local w_u, w_w, lx, lz

    tRamp = timeRamp(t, startRamp[1], startRamp[2])

    FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)    

    t1s = FΓ ⋅ T1s_cs
    t1m2 = t1s ⋅ t1s    
    #t1 = t1s / ((t1s ⋅ t1s).^0.5)

    sΛ = (t1m2.^0.5) / T1m_cs    
    
    vr = UCur_cs*tRamp + waveVel_cs - v
    vn = vr - (vr ⋅ t1s) * t1s / t1m2
    # vn = (v ⋅ t1s) * t1s / t1m2 - v 
    vnm = (vn ⋅ vn).^0.5

    return D_dn * vn * vnm * sΛ 

  end

  function drag_t_ΓX(t, v, u)

    local FΓ, t1s, t1m2, sΛ, vr, vt, vtm, tRamp

    tRamp = timeRamp(t, startRamp[1], startRamp[2])

    FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)

    t1s = FΓ ⋅ T1s_cs
    t1m2 = t1s ⋅ t1s    
    #t1 = t1s / ((t1s ⋅ t1s).^0.5)

    sΛ = (t1m2.^0.5) / T1m_cs
    
    vr = UCur_cs*tRamp + waveVel_cs - v
    vt = (vr ⋅ t1s) * t1s / t1m2
    # vt = -(v ⋅ t1s) * t1s / t1m2    
    vtm = (vt ⋅ vt).^0.5

    return D_dt * vt * vtm * sΛ 

  end
  # ----------------------End----------------------  


  ## Function form added-mass
  # ---------------------Start---------------------
  @unpack C_an, C_at, = params
  D_an = ρw * C_an
  D_at = ρw * C_at   

  printTer("[VAL] D_an = ", D_an)
  printTer("[VAL] D_at = ", D_at)
  printTer()

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

  # Form 0: Simplest
  res0(u, ψu) =          
    ∫( ( (∇(ψu)' ⋅ QTrans_cs) ⊙ (stressK_fnc∘(QTrans_cs, P_cs, ∇(u) )) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ FWeih_cs )*JJ_cs )dΩ 
  
  
  # Form 2: Default
  res(u, ψu) =      
    ∫( ( (∇(ψu)' ⋅ QTrans_cs) ⊙ stressK_fnc(u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ FWeih_cs )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ VectorValue(0.0,1.0) * bedSpring_fnc(u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ drag_n_ΓX(0.0, VectorValue(0.0,0.0), u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ drag_t_ΓX(0.0, VectorValue(0.0,0.0), u) )*JJ_cs )dΩ 


  op_S = FEOperator(res0, US, Ψu)
  # op_S = FEOperator(res, US, Ψu)
  # ----------------------End----------------------


  ## Weak form: Dynamic
  # ---------------------Start---------------------
  
  # Form 0: Simplest
  massD0(t, ∂ₜₜu, ψu) =  
    ∫( ( (ψu ⋅ ∂ₜₜu) * ρcDry )*JJ_cs )dΩ
  # massD(t, u, ∂ₜₜu, v) = massD(t, ∂ₜₜu, v)
  
  resD0(t, u, ψu) =      
    ∫( ( (∇(ψu)' ⋅ QTrans_cs) ⊙ (stressK_fnc∘(QTrans_cs, P_cs, ∇(u) )) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ FWeih_cs )*JJ_cs )dΩ 
    # ∫( (  )*JJ_cs )dΩ      


  # Form 1
  # resD(t, u, ψu) =  
  #   ∫( ( (ψu ⋅ ∂tt(u)) * ρcDry )*JJ_cs )dΩ +
  #   # ∫( ( -ψu ⋅ addedMass_n_ΓX(∂tt(u), u) )*JJ_cs )dΩ +
  #   # ∫( ( -ψu ⋅ addedMass_t_ΓX(∂tt(u), u) )*JJ_cs )dΩ +
  #   ∫( ( -ψu ⋅ VectorValue(0.0,1.0) * bedDamp_fnc(∂t(u), u) )*JJ_cs )dΩ +
  #   ∫( ( (∇(ψu)' ⋅ QTrans_cs) ⊙ stressK_fnc(u) )*JJ_cs )dΩ +
  #   ∫( ( -ψu ⋅ FWeih_cs )*JJ_cs )dΩ +
  #   ∫( ( -ψu ⋅ VectorValue(0.0,1.0) * bedSpring_fnc(u) )*JJ_cs )dΩ +
  #   ∫( ( -ψu ⋅ drag_n_ΓX(t, ∂t(u), u) )*JJ_cs )dΩ +
  #   ∫( ( -ψu ⋅ drag_t_ΓX(t, ∂t(u), u) )*JJ_cs )dΩ 
  #   # ∫( ( -ψu ⋅ drag_ΓX(∂t(u), u) )*JJ_cs )dΩ 
  #   # ∫( (  )*JJ_cs )dΩ   

  # op_D = TransientFEOperator(resD, U, Ψu; order=2)


  # Form 2: Default, No added mass
  massD(t, ∂ₜₜu, ψu) =  
    ∫( ( (ψu ⋅ ∂ₜₜu) * ρcDry )*JJ_cs )dΩ
  # massD(t, u, ∂ₜₜu, v) = massD(t, ∂ₜₜu, v)
  
  resD(t, u, ψu) =      
    # ∫( ( -ψu ⋅ addedMass_n_ΓX(∂tt(u), u) )*JJ_cs )dΩ +
    # ∫( ( -ψu ⋅ addedMass_t_ΓX(∂tt(u), u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ VectorValue(0.0,1.0) * bedDamp_fnc(∂t(u), u) )*JJ_cs )dΩ +
    ∫( ( (∇(ψu)' ⋅ QTrans_cs) ⊙ stressK_fnc(u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ FWeih_cs )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ VectorValue(0.0,1.0) * bedSpring_fnc(u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ drag_n_ΓX(t, ∂t(u), u) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ drag_t_ΓX(t, ∂t(u), u) )*JJ_cs )dΩ 
    # ∫( ( -ψu ⋅ drag_ΓX(∂t(u), u) )*JJ_cs )dΩ 
    # ∫( (  )*JJ_cs )dΩ     


  # Define operator
  op_D = TransientSemilinearFEOperator(massD0, resD0, U, Ψu; 
    order=2, constant_mass=true)
  # op_D = TransientSemilinearFEOperator(massD, resD, U, Ψu; 
  #   order=2, constant_mass=true)

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

  printTer("[RES] Length Catenary solution = ", calcLen(J) )  
  printTer("[RES] Length Static Solution = ", calcLen(JNew(uh_S)) )  
  printTer()

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

  # Implicit solver
  ode_solver = GeneralizedAlpha2(nls, simΔt, 0.0)
  # GenAlpha 1.0 Midpoint: Diverges quickly
  # GenAlpha 0.0 Fully implicit: Stable
  
  solnht = solve(ode_solver, op_D, t0, simT, (U0,U0t)) 
  # solnht = solve(ode_solver, op_D, t0, simT, (U0,U0t,U0tt)) 
  # ----------------------End----------------------


  ## Save quantities
  # ---------------------Start---------------------  
  rPrb = Point.(0.0:L/5:L)
  nPrb = length(rPrb)

  daFile1 = open( pltName*"data1.dat", "w" )
  
  # ----------------------End----------------------


  pvd = paraview_collection(pltName*"tSol", append=false)
  if(outFreeSurface)
    if(isdir(pltName*"fs/"))
      rm(pltName*"fs/", recursive = true)
    end
    mkpath(pltName*"fs/")
    pvd_fs = paraview_collection(pltName*"fs/fs_tSol", append=false)
  end


  ## Save initial solution
  # ---------------------Start---------------------  
  uh = U0    
  @printf("Time : %10.3f s \t Counter : %5i \n", t0, 0)          
  tprt = @sprintf("%d",floor(Int64,t0*1000))                        

  xNew = X + uh

  pvd[t0] = createvtk(Ω,    
    pltName*"tSol_$tprt"*".vtu",
    cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh, 
      "ETang"=>ETang(uh), "sigma"=>stressσ(uh),
      "spr"=>spng(uh) ])
      # "drag_ΓX" => drag_ΓX_intp(uh,uh) ])

  if(outFreeSurface)
    pvd_fs[t0] =createvtk(Ω_fs,    
      pltName*"fs/fs_tSol_$tprt"*".vtu",
      cellfields=["X"=>X_fs, "eta"=> getEta_fs(t0)])
  end
  
  # ----------------------End----------------------


  ## Execute
  # ---------------------Start---------------------

  # Interpolation: save cache
  xNew = X + uh
  save_cache1, save_cache2 = Gridap.Arrays.return_cache(xNew, rPrb)
  save_f_cache = save_cache2[2]
  

  execTime = zeros(Float64, 1, 10)
  execTime[1] = time()  # tick()
  execTime[3] = time()   
  tick()
  cnt=0
  
  # Update waveVel(tn) before solving t(n+1)
  update_state!( (a,b) -> (true, b), waveVel_cs, 
    getWaveVel_cf(t0, xNew) ) 
  
  # for (t, uh) in solnht                       
  next = iterate(solnht)            
  while next !== nothing

    (iSol, iState) = next
    (t, uh) = iSol      
    iNLCache = iState[2][5][1][4]
    # @show propertynames(iNLCache.result)

    cnt = cnt+1          
    @printf("Time : %10.3f s \t Counter : %5i \n", t, cnt)          
    @printf("Conv : %10s \t Iter    : %5i \n",
      iNLCache.result.x_converged, iNLCache.result.iterations)
    tprt = @sprintf("%d",floor(Int64,t*1000000))

    # Interpolation: easiest method
    # xNew = X + uh
    # σT = getPrincipalStress(uh, rPrb)      
    # xNewPrb = xNew.(rPrb)
    
    # Interpolation: using return_cache()
    # xNew = X + uh
    # cache_xNew = Gridap.Arrays.return_cache(xNew, rPrb)
    # xNewPrb = evaluate!(cache_xNew, xNew, rPrb)

    # Interpolation: open return_cache(), re-use sub_cache
    xNew = X + uh   
    cell_f = get_array(xNew)
    cell_f_cache = array_cache(cell_f)    
    cache2 = cell_f_cache, save_f_cache, cell_f, xNew
    cache_xNew = (save_cache1, cache2)
    xNewPrb = evaluate!(cache_xNew, xNew, rPrb)
    

    sT_uh = stressσ_fnc∘(QTrans, P, J, ∇(uh) )
    cell_f = get_array(sT_uh)
    cell_f_cache = array_cache(cell_f)
    cache2 = cell_f_cache, save_f_cache, cell_f, sT_uh
    cache_sT = (save_cache1, cache2)
    sTPrb = evaluate!(cache_sT, sT_uh, rPrb)
    σT = getPrincipalStress2.(sTPrb)      


    @printf(daFile1, "%15.3f",t)
    @printf(daFile1, ", %2i, %5i", 
      iNLCache.result.x_converged, iNLCache.result.iterations)
    # [print(daFile1, string(val)*", \t") for val in lDa]
    [@printf(daFile1, ", %15.3f, %15.3f, %15.3f, %15.3f", 
      rPrb[i][1], xNewPrb[i][1], xNewPrb[i][2], σT[i]*A_str)
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

      if(outFreeSurface)
        pvd_fs[t] =createvtk(Ω_fs,    
          pltName*"fs/fs_tSol_$tprt"*".vtu",
          cellfields=["X"=>X_fs, "eta"=> getEta_fs(t)])
      end
    
    end


    execTime[4] = time()  
    tock()
    @printf(outFile1, 
      "%5i, %10.3f, %10.3f, %5i, %2i \n", 
      cnt, t, execTime[4]-execTime[3], 
      iNLCache.result.iterations, iNLCache.result.x_converged)
    println("-x-x-x-")
    println()
    execTime[3] = time()  
    tick()

    # # Update waveVel(tn) before solving t(n+1)
    # update_state!( (a,b) -> (true, b), waveVel_cs, 
    #   getWaveVel_cf(t, xNew) ) 
    
    next = iterate(solnht, iState)
  end  
  execTime[2] = time()  
  tock()
  @printf(outFile0, 
    "\n[TIM] Total Time: \t %5i \t %.3f \n", 
    round(simT/simΔt), execTime[2]-execTime[1])

  vtk_save(pvd)
  if(outFreeSurface)
    vtk_save(pvd_fs)
  end

  close(outFile0)
  close(outFile1)
  close(daFile1)
  # ----------------------End----------------------

    
end




"""
Aux functions
=============

"""
# ---------------------Start---------------------

function printTerAndFile(str::String, 
    val::Union{AbstractArray,Tuple,Real}, outFile::IOStream)
  
  println(str, val)
  println(outFile, str, val)
  # @printf("%s %15.6f\n", str, val)
  # @printf(outFile,"%s %15.6f\n", str, val)
end

function printTerAndFile(str::String, outFile::IOStream)
  println(str)
  println(outFile, str)
end


function getInputSpec(params)

  @unpack Hs, Tp, h0, nω, seed, ωc = params

  if(ωc < 0)
    ω, S, A = jonswap(Hs, Tp,
      plotflag=false, nω = nω)
  else
    ω, S, A = jonswap(Hs, Tp,
      plotflag=false, nω = nω, ωc = ωc)
  end

  k = dispersionRelAng.(h0, ω; msg=false)
  α = randomPhase(ω, seed = seed)

  sp = SpecStruct( h0, ω, S, A, k, α; Hs = Hs, Tp = Tp )
  return sp
end


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