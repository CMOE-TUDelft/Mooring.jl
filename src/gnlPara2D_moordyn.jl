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
using Plots

using Mooring.BedSpring
using Mooring.Drag
using Mooring.GnlCommon
using Mooring.StressLinear

include(srcdir("subroutines","testParams.jl"))



"""
main()
======

"""
function main(params)

  @unpack resDir = params
  pltName = resDir*"/gnl_"


  ## Terminal Output
  # ---------------------Start---------------------  
  outFile0 = open( pltName*"aout.dat", "w" )
  outFile1 = open( pltName*"arunTime.dat", "w" )
  
  printTer(a::String,b) = printTerAndFile(a,b,outFile0)
  printTer(a::String) = printTerAndFile(a,outFile0)
  printTer() = printTerAndFile("",outFile0)
  # ----------------------End----------------------  


  @unpack ρw = params #Kg/m3 density of water    

  # Line properties
  @unpack ϵ0, xz_fl = params  
  seg = StressLinear.Segment( params )
  printTer("[VAL] Length = ", seg.L)
  printTer("[VAL] Given FL xz = ", xz_fl)
  printTer()  

  @unpack bedObj = params  
  bedObj = BedSpring.setStillWei!(bedObj, seg.ρcSub*g)    
  printTer("[VAL] Bed StillWei = ", bedObj.stillWei)
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
  @unpack nx, order = params
  domain = (0, seg.L)
  partition = (nx)      
  model = CartesianDiscreteModel(domain, partition)

  # Labelling 
  labels_Ω = get_face_labeling(model)
  add_tag_from_tags!(labels_Ω,"anchor",[1]) 
  add_tag_from_tags!(labels_Ω,"fairLead",[2]) 
  writevtk(model, pltName*"model")

  printTer("[VAL] nx = ", nx)
  printTer("[VAL] order = ", order)
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
  model_fs = CartesianDiscreteModel((0, 1.2*seg.L), (nx))
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


  # ## Reference config Catenary
  # # ---------------------Start---------------------    
  # interpX, interpZ = setInitXZ(initCSV)

  # X(r) = VectorValue( interpX(r[1]/L), interpZ(r[1]/L) )

  # writevtk(Ω, pltName*"referenceDomain",
  #   cellfields=["X"=>X])
  # # ----------------------End----------------------  


  ## Reference config parabola
  # ---------------------Start---------------------        
  X, printMsg = getParabola(xz_fl[1],xz_fl[2], seg.L)  
  printTer(printMsg)

  writevtk(Ω, pltName*"referenceDomain",
    cellfields=["X"=>X])
  # ----------------------End----------------------  


  ## Define Test Fnc
  # ---------------------Start---------------------
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
  gFairLead_S(x) = VectorValue(0.0, ϵ0*seg.L)  

  US = TrialFESpace(Ψu, [gAnch_S, gFairLead_S])
  # ----------------------End----------------------


  ## Define Trial Fnc Dynamic
  # ---------------------Start---------------------  
  # Dirichlet BC
  gAnch(x, t::Real) = VectorValue(0.0, 0.0)
  gAnch(t::Real) = x -> gAnch(x, t)

  Xh_fl = X(Point(seg.L))
  printTer("[VAL] Xh_fl = ", (Xh_fl[1], Xh_fl[2]))
  printTer("[VAL] Xh_fl = ", (Xh_fl[1], Xh_fl[2]-h0))
  printTer()  
  
  @unpack ffm_η, ffm_ω = params
  printTer("[VAL] ffm_η ", ffm_η)
  printTer("[VAL] ffm_ω ", ffm_ω)
  function getFairLeadEnd(x,t)    
    # ffm_η = 1.5 #m2
    # ffm_ω = 1 #rad/s

    tRamp = timeRamp(t, startRamp[1], startRamp[2])    

    return VectorValue( 
      tRamp*ffm_η*sin(ffm_ω*t),
      ϵ0*seg.L )
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
  FWeih = interpolate_everywhere(VectorValue(0.0, -seg.ρcSub*g), Ψu)
  
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
  stressK(u) = 2*seg.μm * (FΓ(u) ⋅ ETang(u))
  stressS(u) = 2*seg.μm * ETang(u)
  stressσ(u) = ( FΓ(u) ⋅ stressS(u) ⋅ FΓ(u)' ) / sΛ(u) 
    

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
  


  ## Function form drag
  # ---------------------Start---------------------
  @unpack C_dn, d_dn, C_dt, d_dt = params  
  D_dn = 0.5 * ρw * C_dn * d_dn / seg.A #kg/m4
  D_dt = 0.5 * ρw * C_dt * π * d_dt / seg.A #kg/m4
  
  printTer("[VAL] D_dn = ", D_dn)
  printTer("[VAL] D_dt = ", D_dt)
  printTer()
  
  
  function drag_ΓX(QTr, T1s, T1m, ∇u, v) #No wave and current

    local FΓ, t1s, t1m2, vn, vnm, sΛ, vt, vtm

    FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)

    t1s = FΓ ⋅ T1s
    t1m2 = t1s ⋅ t1s    
    #t1 = t1s / ((t1s ⋅ t1s).^0.5)

    sΛ = (t1m2.^0.5) / T1m
    
    vt = -(v ⋅ t1s) * t1s / t1m2
    vtm = (vt ⋅ vt).^0.5
    vn = -v - vt
    vnm = (vn ⋅ vn).^0.5    

    return (D_dn * vn * vnm + D_dt * vt * vtm) * sΛ 

  end    
  # ----------------------End----------------------



  ## Tuples for cellstates
  csTup1 = (QTrans_cs, T1s_cs, T1m_cs)



  ## Parsing functions
  # ---------------------Start---------------------
  bedSpring_fnc(X, QTr, T1s, T1m, u, ∇u, v) = 
    BedSpring.forceFnc(bedObj, X, QTr, T1s, T1m, u, ∇u, v)

  stressK_fnc(QTr, P, ∇u) = 
    StressLinear.stressK_fnc(seg, QTr, P, ∇u)
  
  stressσ_fnc(QTr, P, J, ∇u ) = 
    StressLinear.stressσ_fnc(seg, QTr, P, J, ∇u )
  
  ETang_fnc(QTr, P, J, ∇u ) = 
    StressLinear.ETang_fnc( QTr, P, J, ∇u )
  # ----------------------End----------------------



  ## Weak form: Static 
  # ---------------------Start---------------------

  # Form 0: Simplest
  res0(u, ψu) =          
    ∫( ( (∇(ψu)' ⋅ QTrans_cs) ⊙ (stressK_fnc∘(QTrans_cs, P_cs, ∇(u) )) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ FWeih_cs )*JJ_cs )dΩ + 
    ∫( ( -ψu ⋅ VectorValue(0.0,1.0) * 
      (bedSpring_fnc∘(Xh_cs, csTup1...,u, ∇(u), 0.0*u)) )*JJ_cs )dΩ 


  op_S = FEOperator(res0, US, Ψu)
  # op_S = FEOperator(res, US, Ψu)
  # ----------------------End----------------------


  ## Weak form: Dynamic
  # ---------------------Start---------------------    

  # Form 1: Self drag only
  massD1(t, ∂ₜₜu, ψu) =  
    ∫( ( (ψu ⋅ ∂ₜₜu) * seg.ρcDry )*JJ_cs )dΩ
  # massD(t, u, ∂ₜₜu, v) = massD(t, ∂ₜₜu, v)
  
  resD1(t, u, ψu) =      
    ∫( ( (∇(ψu)' ⋅ QTrans_cs) ⊙ (stressK_fnc∘(QTrans_cs, P_cs, ∇(u))) )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ FWeih_cs )*JJ_cs )dΩ +
    ∫( ( -ψu ⋅ VectorValue(0.0,1.0) * 
      (bedSpring_fnc∘(Xh_cs, csTup1..., u, ∇(u), ∂t(u))) )*JJ_cs )dΩ +    
    ∫( ( -ψu ⋅ (drag_ΓX∘(csTup1..., ∇(u), ∂t(u))) )*JJ_cs )dΩ 
    # ∫( (  )*JJ_cs )dΩ        

    
  # Define operator
  # op_D = TransientSemilinearFEOperator(massD0, resD0, U, Ψu; 
  #   order=2, constant_mass=true)
  op_D = TransientSemilinearFEOperator(massD1, resD1, U, Ψu; 
    order=2, constant_mass=true)
  # ----------------------End----------------------


  ## Static Solution
  # ---------------------Start---------------------
  tick()
  fx_U(r) = 0.01*sin(π*r[1])
  fy_U(r) = 0.0*sin(2*π*r[1])
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
  tock()
  printTer()

  xNew = Xh + uh_S

  writevtk(Ω, pltName*"staticRes",
    cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh_S, 
      "ETang"=>ETang_fnc( QTrans, P, J, ∇(uh_S) ), 
      "sigma"=>stressσ_fnc( QTrans, P, J, ∇(uh_S) ) ])
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
  # GenAlpha is always stable
  # GenAlpha 1.0 Midpoint 
  #   No dissipation case: Can Diverge due to high freq
  # GenAlpha 0.0 Fully implicit
  #   Asymptotic annhilition: Highly dissipative
  #   T < 10*Δt is dissipated
  # GenAlpha 0.4 
  #   Used in OrcaFlex implicit
  
  solnht = solve(ode_solver, op_D, t0, simT, (U0,U0t)) 
  # solnht = solve(ode_solver, op_D, t0, simT, (U0,U0t,U0tt)) 
  # ----------------------End----------------------


  ## Save quantities
  # ---------------------Start---------------------  
  rPrbMat = [0:seg.L/20:seg.L;]
  # rPrb = Point.(0.0:L/10:L)
  rPrb = Point.(rPrbMat)
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
      "ETang"=>ETang_fnc( QTrans, P, J, ∇(uh) ), 
      "sigma"=>stressσ_fnc( QTrans, P, J, ∇(uh) ),
      "gradU"=>∇(uh)])      

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
  save_f_cache2 = save_cache2[2]
  

  execTime = zeros(Float64, 1, 10)
  execTime[1] = time()  # tick()
  execTime[3] = time()   
  tick()
  cnt=0
  
  
  # for (t, uh) in solnht                       
  next = iterate(solnht)            
  while next !== nothing

    (iSol, iState) = next
    (t, uh) = iSol      
    iNLCache = iState[2][5][1][4]
    # @show propertynames(iNLCache.result)

    cnt = cnt+1          
    @printf("Progress : %10.3f %% \n", t/simT*100)          
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
    cache2 = assemble_cache(xNew, save_f_cache2)
    xNewPrb = evaluate!((save_cache1, cache2), xNew, rPrb)
    

    sT_uh = stressσ_fnc∘(QTrans, P, J, ∇(uh) )
    cache2 = assemble_cache(sT_uh, save_f_cache2)
    sTPrb = evaluate!((save_cache1, cache2), sT_uh, rPrb)
    sTPrb_princ = getPrincipalStress2.(sTPrb)


    ETang_uh = ETang_fnc∘(QTrans, P, J, ∇(uh) )
    cache2 = assemble_cache(ETang_uh, save_f_cache2)
    ETangPrb = evaluate!((save_cache1, cache2), ETang_uh, rPrb)
    ETangPrb_princ = getPrincipalStress2.(ETangPrb)


    # gradU = ∇(uh)
    # cache2 = assemble_cache(gradU, save_f_cache2)
    # gradUPrb = evaluate!((save_cache1, cache2), gradU, rPrb)
    # @show gradUPrb[1]
    

    @printf(daFile1, "%15.6f",t)
    @printf(daFile1, ", %2i, %5i", 
      iNLCache.result.x_converged, iNLCache.result.iterations)
    # [print(daFile1, string(val)*", \t") for val in lDa]
    [@printf(daFile1, ", %15.6f, %15.6f, %15.6f, %20.10e, %20.10e", 
      rPrb[i][1], xNewPrb[i][1], xNewPrb[i][2], 
      ETangPrb_princ[i], sTPrb_princ[i])
      for i in 1:nPrb]
    @printf(daFile1, "\n")

          
    if(cnt%outMod == 0)               

      println("Paraview output")        

      pvd[t] = createvtk(Ω,    
        pltName*"tSol_$tprt"*".vtu",
        cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh, 
          "ETang"=>ETang_fnc( QTrans, P, J, ∇(uh) ), 
          "sigma"=>stressσ_fnc( QTrans, P, J, ∇(uh) ),
          "gradU"=>∇(uh)])

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



  daF = CSV.read(pltName*"data1.dat", DataFrame, header=false)
  daF = Matrix(daF)
  fT = 2*pi/ffm_ω


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

  nVals = 5

  prbi = nPrb-1
  plt1 = plot()
  plot!( plt1, daF[:,1]/fT, daF[:,3+(prbi-1)*nVals+2], 
    label = nothing, linewidth=3 )
  plot!( plt1, title = "X Coord (FL-1)", ylabel = "(m)" )
  setGridlines(plt1)

  prbi = nPrb-1
  plt2 = plot()
  plot!( plt2, daF[:,1]/fT, daF[:,3+(prbi-1)*nVals+3], 
    label = nothing, linewidth=3 )
  plot!( plt2, title = "Z Coord (FL-1)", ylabel = "(m)" )
  setGridlines(plt2)

  prbi = 1
  plt3 = plot()
  plot!( plt3, daF[:,1]/fT, daF[:,3+(prbi-1)*nVals+4], 
    label = nothing, linewidth=3 )
  plot!( plt3, title = "ETang at Anchor", ylabel = "ETang" )
  setGridlines(plt3)

  prbi = nPrb
  plt4 = plot()
  plot!( plt4, daF[:,1]/fT, daF[:,3+(prbi-1)*nVals+4], 
    label = nothing, linewidth=3 )
  plot!( plt4, title = "ETang at FairLead", ylabel = "ETang" )
  setGridlines(plt4)


  prbi = nPrb-1
  plt5 = plot()
  plot!( plt5, daF[:,1]/fT, daF[:,3+(prbi-1)*nVals+4], 
    label = nothing, linewidth=3 )
  plot!( plt5, title = "ETang at FairLead - 1", ylabel = "ETang" )
  setGridlines(plt5)

  savefig(plt1, pltName*"plot_posX.png")
  savefig(plt2, pltName*"plot_posZ.png")
  savefig(plt3, pltName*"plot_ETang_Anchor.png")
  savefig(plt4, pltName*"plot_ETang_FairLead.png")
  savefig(plt5, pltName*"plot_ETang_FairLeadNear.png")

    
end


end