
@with_kw struct gnl_simple_params
  η₀::Float64 = 1.0 # wave amplitude
  ω::Float64 = 0.1 # wave frequency
  E::Float64 = 1e4 # rope Young's modulus
  ρₘ::Float64 = 1 # rope density
  L::Float64 = 1 # rope length
  A::Float64 = 0.01 # rope cross-sectional area
  nx::Int = 100 # number of elements
  verbose::Bool = false # print debug info
  Δt::Float64 = 0.1 # time step
  T::Float64 = 0.2 # final time
  outΔt::Float64 = 0.1 # output time step
  testname::String = "tmp" # test name
  sampling_points::Vector{Float64} = [0.5] # sampling point (relative to the length)
end

function main_dynamic(params)

  # Data filename
  @unpack testname = params
  filename = datadir(testname)
  csvfile = open(filename*".csv", "w")

  # Material properties
  @unpack E,ρₘ,L,A = params
  μₘ = 0.5*E
  EA = E*A
  m = ρₘ * A
  L

  # Parameter Domain
  @unpack nx = params
  domain = (0, L)
  partition = (nx)
  model = CartesianDiscreteModel(domain, partition)

  # Labelling
  @unpack verbose = params
  labels_Ω = get_face_labeling(model)
  add_tag_from_tags!(labels_Ω,"leftEdge",[1])
  add_tag_from_tags!(labels_Ω,"rightEdge",[2])
  verbose && writevtk(model, filename*"_model")

  # Triangulations
  Ω = Interior(model)
  Γ = Boundary(model)

  ## Reference config Fries paper 5.2
  # ---------------------Start---------------------
  fx_X(r) = r[1]
  fy_X(r) = 1/2*(1-r[1]) - 1/7*sin(π*r[1])
  X(r) = VectorValue(fx_X(r), fy_X(r))
  verbose && writevtk(Ω, filename*"_referenceDomain",
    cellfields=["X"=>X])
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
  nΓ = get_normal_vector(Γ)

  ## Geometric quantities
  # ---------------------Start---------------------
  Xh = interpolate_everywhere(X, Ψu)

  J = ∇(Xh)'
  LOrig = sum(∫( sqrt∘(J ⊙ J) )dΩ)
  verbose && println("Reference Length = ", LOrig )
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

  # External loads
  @unpack η₀,ω = params
  g=9.81
  FBodyh(t) = VectorValue(η₀*sin(ω*t), -ρₘ*g)#interpolate_everywhere(VectorValue(η₀*sin(ω*t), -ρₘ*g), Ψu)
  bedK = ρₘ*g*2
  bedRamp = 1e3
  spng(u) = 0.5+0.5*(tanh∘( VectorValue(0.0,-bedRamp) ⋅ (Xh+u)))

  ## Weak form
  # ---------------------Start---------------------
  mass(t, u, ψu) =
    ∫( ( (ψu ⋅ ∂tt(u)) * ρₘ ) * ( (J ⊙ J).^0.5 ) )dΩ

  res(t,u, ψu) =
    ∫(
      (
        ∇X_Dir(ψu) ⊙ stressK(u) +
        - ( ψu ⋅ FBodyh(t) ) +
        - ( ψu ⋅ VectorValue(0.0,bedK) * spng(u) )
      )*((J ⊙ J).^0.5)
    )dΩ

  resD(t, u, ψu) = mass(t,u,ψu) + res(t,u,ψu)

  op_AD = TransientFEOperator(resD, U, Ψu; order=2)
  # ----------------------End----------------------

  ## Initial solution
  # ---------------------Start--------------------
  t0 = 0.0
  da = wload(datadir(filename*"_sol0.jld2"))
  uf = FEFunction( Ψu, da["uh"] )
  U0 = interpolate_everywhere(uf, U(t0))
  U0t = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))
  U0tt = interpolate_everywhere(VectorValue(0.0, 0.0), U(t0))
  # ----------------------End----------------------

  ## Solver
  # ---------------------Start---------------------
  @unpack Δt, T = params
  nls = NLSolver(show_trace=true,
    method=:newton, linesearch=Static(),
    iterations=100, ftol = 1e-8)
  ode_solver = GeneralizedAlpha(nls, Δt, 0.0)
  solnht = solve(ode_solver, op_AD, (U0,U0t,U0tt), t0, T)
  # ----------------------End----------------------

  # Execute
  @unpack outΔt, sampling_points = params
  outMod = floor(Int64,outΔt/Δt)
  points = [Point(L*p,) for p in sampling_points]
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

      Sₕ = stressS(uh)
      Sₕ_norms = [norm(Sₕ(point)) for point in points]
      vals = [t, Sₕ_norms...]
      CSV.write(csvfile, Tables.table(vals'), append=true)

      if (verbose)
        pvd[t] = createvtk(Ω,
          filename*"_tSol_$tprt"*".vtu",
          cellfields=["XOrig"=>X, "XNew"=>xNew, "uh"=>uh,
          "ETang"=>ETang(uh), "sigma"=>stressσ(uh),
          "spr"=>spng(uh) ])
      end
    end
  end

  close(csvfile)

end

function main_static(params)

  # Data filename
  @unpack testname = params
  filename = datadir(testname)

  # Material properties
  @unpack E,ρₘ,L,A = params
  μₘ = 0.5*E
  EA = E*A
  m = ρₘ * A
  L

  # Parameter Domain
  @unpack nx = params
  domain = (0, L)
  partition = (nx)
  model = CartesianDiscreteModel(domain, partition)

  # Labelling
  @unpack verbose = params
  labels_Ω = get_face_labeling(model)
  add_tag_from_tags!(labels_Ω,"leftEdge",[1])
  add_tag_from_tags!(labels_Ω,"rightEdge",[2])
  verbose && writevtk(model, filename*"_model")

  # Triangulations
  Ω = Interior(model)
  Γ = Boundary(model)

  ## Reference config Fries paper 5.2
  # ---------------------Start---------------------
  fx_X(r) = r[1]
  fy_X(r) = 1/2*(1-r[1]) - 1/7*sin(π*r[1])
  X(r) = VectorValue(fx_X(r), fy_X(r))
  verbose && writevtk(Ω, filename*"_referenceDomain",
    cellfields=["X"=>X])
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
  nΓ = get_normal_vector(Γ)


  ## Geometric quantities
  # ---------------------Start---------------------
  Xh = interpolate_everywhere(X, Ψu)

  J = ∇(Xh)'
  LOrig = sum(∫( sqrt∘(J ⊙ J) )dΩ)
  verbose && println("Reference Length = ", LOrig )
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

  # External loads
  @unpack η₀,ω = params
  g = 9.81
  FBodyh = VectorValue(0.0, -ρₘ*g)#interpolate_everywhere(VectorValue(η₀*sin(ω*t), -ρₘ*g), Ψu)
  bedK = ρₘ*g*2
  bedRamp = 1e3
  spng(u) = 0.5+0.5*(tanh∘( VectorValue(0.0,-bedRamp) ⋅ (Xh+u)))

  ## Weak form
  # ---------------------Start---------------------
  res(u, ψu) =
  ∫( (
      ∇X_Dir(ψu) ⊙ stressK(u) +
      - ( ψu ⋅ FBodyh ) +
      - ( ψu ⋅ VectorValue(0.0,bedK) * spng(u) )
    )*((J ⊙ J).^0.5) )dΩ
    op = FEOperator(res, U, Ψu)
  # ----------------------End----------------------

  ## Initial solution
  # ---------------------Start--------------------
  fx_U(r) = 0.0#sin(π*r[1])
  fy_U(r) = -0.0000001*sin(π*r[1])
  Ua(r) = VectorValue(fx_U(r), fy_U(r))
  # Ua(r) = VectorValue(0.0, -0.00001)
  U0 = interpolate_everywhere(Ua, U)
  # ----------------------End----------------------

  ## Solver
  # ---------------------Start---------------------
  nls = NLSolver(show_trace=true,
    method=:newton, linesearch=Static(),
    iterations=100, ftol = 1e-8)
  (uh, cache) = solve!(U0, nls, op)
  # ----------------------End----------------------

  data = Dict(
    "uh" => get_free_dof_values(uh)
  )

  wsave(filename*"_sol0.jld2", data)

end
