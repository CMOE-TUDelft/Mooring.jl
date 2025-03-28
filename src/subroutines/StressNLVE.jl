module StressNLVE

using Revise
using Gridap
using ForwardDiff
using Parameters
using Roots: find_zero

using Mooring.StressLinear: Segment

ITER_SOLVE = 1
ITER_DIFF = 1

"""
Custom Structs
=============

"""
# ---------------------Start---------------------
@with_kw struct Schapery

  D0::Real
  N::Real
  Dn::Vector{Real}
  λn::Vector{Real}
  g0::Tuple{Bool,Real,Vector{Real}} #Vector{Float64}  
  g1::Tuple{Bool,Real,Vector{Real}} #Vector{Float64}  
  g2::Tuple{Bool,Real,Vector{Real}} #Vector{Float64}  
  # Tuple{Bool,Real,Vector{Real}}
  # Bool: Turn on limiter?
  # Real: Limiter value
  # Vector{Real}: Coefficients

  # Calculated values
  ΣDn::Real

end


function Schapery(
  placeholder; # Ask Oriol # TODO
  D0::Real,
  Dn::Vector{<:Real}=[0.0],
  λn::Vector{<:Real}=[1.0],
  g0::Union{ Tuple{Real, Vector{<:Real}}, Vector{<:Real} } = [1.0],
  g1::Union{ Tuple{Real, Vector{<:Real}}, Vector{<:Real} } = [1.0],
  g2::Union{ Tuple{Real, Vector{<:Real}}, Vector{<:Real} } = [1.0]
)
  # Determine the type of g0, g1, g2 and normalize them to Tuple form
  g0_tuple = isa(g0, Tuple) ? (true, g0[1], g0[2]) : (false, 0.0, g0)
  g1_tuple = isa(g1, Tuple) ? (true, g1[1], g1[2]) : (false, 0.0, g1)
  g2_tuple = isa(g2, Tuple) ? (true, g2[1], g2[2]) : (false, 0.0, g2)

  nDn = length(Dn)
  nλn = length(λn)

  N = min(nDn, nλn)
  ΣDn = sum(Dn[1:N])

  Schapery(D0, N,
    Dn[1:N], λn[1:N],
    g0_tuple, g1_tuple, g2_tuple,
    ΣDn)
end



mutable struct SchaperyData
  qt0
  qt1
  pETang_t0
  pETang_t1
  pS_t0
  pS_t1

  function SchaperyData()
    new()
  end
end
# ----------------------End----------------------



"""
Functions
=============

"""
# ---------------------Start---------------------

## Function to rotate to get principal strain / stress
function getStrRotMatrix(σ)
  # Ref
  # https://www.continuummechanics.org/principalstress.html
  
  θp = atan( 2*σ[2]/(σ[1]-σ[4]) ) / 2
  return TensorValue(cos(θp), -sin(θp), sin(θp), cos(θp))
end


## General polynomial evaluator
function poly_eval(
  σ, 
  coeffsAll::Tuple{Bool,Real,Vector{Real}}
)
  
  limOn, lim, coeffs = coeffsAll

  if limOn      
    if σ > lim
      return sum(c * σ^(i - 1) for (i, c) in enumerate(coeffs))
    else
      return 1.0  # Linear viscoelastic behavior below threshold
    end
  end

  return sum(c * σ^(i - 1) for (i, c) in enumerate(coeffs))
end
# ----------------------End----------------------



"""
Stress-strain functions
=============

"""
# ---------------------Start---------------------
function stressK_NLVE(sch::Schapery, Δt,
  QTr, P, ∇u, 
  schDa1_ϵt0, schDa1_qt0, schDa1_pS_t0,
  schDa2_ϵt0, schDa2_qt0, schDa2_pS_t0)
    
	local FΓ, EDir, ETang
  local pS_tk1, pS_tk2, err1, err2, pS_tguess
  
	
	FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
	# FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)
	EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )
	ETang = P ⋅ EDir ⋅ P

	# return 2*seg.μm * (FΓ ⋅ ETang)  

  rotM = getStrRotMatrix( ETang )
  pETang = rotM ⋅ (ETang ⋅ transpose(rotM))  

  pS_tguess = schDa1_pS_t0
  pS_tk1, err1 = get_stressNLVE(
    sch, Δt,
    schDa1_ϵt0, schDa1_qt0.data, schDa1_pS_t0,
    pETang[1], pS_tguess )

  pS_tguess = schDa2_pS_t0
  pS_tk2, err2 = get_stressNLVE(
    sch, Δt,
    schDa2_ϵt0, schDa2_qt0.data, schDa2_pS_t0,
    pETang[4], pS_tguess )

  pStr = TensorValue( 
    pS_tk1,
    0.0, 0.0, 
    pS_tk2 )

  # pStr = TensorValue( 
  #     linStressStrain(seg, pETang[1]),
  #     0.0, 0.0, 
  #     0.0 )
  #     # linStressStrain(seg, pETang[4]) )

  S = ( transpose(rotM) ⋅ pStr ) ⋅ rotM

  return FΓ ⋅ S
end


function get_stressNLVE(
  sch, Δt,
  ϵt0, qt0, pS_t0,
  # ϵt1, 
  ϵt1::ForwardDiff.Dual, 
  pS_guess )

  local pS_tk1, err1
  
  # Solving by fixed point iteration
  # TODO: number of iterations 
  pS_tk1 = pS_guess
  for i = 1:ITER_DIFF
    pS_tk1, err1 = StressNLVE.σPredicted( 
      sch, 
      ϵt0, Δt, qt0, pS_t0,
      ϵt1, Δt, pS_tk1 )
  end

  return pS_tk1, err1
end


function get_stressNLVE(
  sch, Δt,
  ϵt0, qt0, pS_t0,
  ϵt1::Float64, 
  pS_guess )

  local pS_tk1, err1
  
  # Solving by fixed point iteration
  # TODO: number of iterations 
  pS_tk1 = pS_guess
  for i = 1:ITER_SOLVE
    pS_tk1, err1 = StressNLVE.σPredicted( 
      sch, 
      ϵt0, Δt, qt0, pS_t0,
      ϵt1, Δt, pS_tk1 )
  end

  # @show "fixed-point iteration"

  return pS_tk1, err1
end



# function get_stressNLVE(
#   sch, Δt,
#   ϵt0, qt0, pS_t0,
#   ϵt1::Float64, pS_guess )

#   local pS_tk1, err1

#   err(σi) = StressNLVE.Residual_σPredicted( 
#     sch, 
#     ϵt0, Δt, qt0, pS_t0,
#     ϵt1, Δt, σi )
  
#   pS_tk1 = find_zero(err, pS_guess)    
#   # @show "find_zero"

#   return pS_tk1, 0.0
# end


function get_rotM(QTr, P, J, ∇u)
    
	local FΓ, EDir, ETang
	
	FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
	EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )

	ETang = P ⋅ EDir ⋅ P

  rotM = getStrRotMatrix( ETang )
  
  return rotM
end


function update_pETang(QTr, P, J, ∇u, index)
    
	local FΓ, EDir, ETang
	
	FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
	# FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)
	EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )

	ETang = P ⋅ EDir ⋅ P

  rotM = getStrRotMatrix( ETang )
  pETang = rotM ⋅ (ETang ⋅ transpose(rotM))  

  return pETang[index]
end


# function update_pS(sch::Schapery, Δt,
#   QTr, P, ∇u, 
#   schDa1_ϵt0, schDa1_qt0, schDa1_pS_t0)
    
# 	local FΓ, EDir, ETang
#   local pS_tk1, err1
	
# 	FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
# 	# FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)
# 	EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )
# 	ETang = P ⋅ EDir ⋅ P

# 	# return 2*seg.μm * (FΓ ⋅ ETang)  

#   rotM = getStrRotMatrix( ETang )
#   pETang = rotM ⋅ (ETang ⋅ transpose(rotM))  

#   pS_tk1 = schDa1_pS_t0
#   for i = 1:1
#     pS_tk1, err1 = StressNLVE.σPredicted( 
#       sch, 
#       schDa1_ϵt0, Δt, schDa1_qt0.data, schDa1_pS_t0,
#       pETang[1], Δt, pS_tk1 )
#   end

#   return pS_tk1    
# end


function update_pS(sch::Schapery, Δt,
  schDa1_ϵt0, schDa1_qt0, schDa1_pS_t0,
  schDa1_ϵt1)
    
  local pS_tk1, pS_tguess, err1

  pS_tguess = schDa1_pS_t0
  pS_tk1, err1 = get_stressNLVE(
    sch, Δt,
    schDa1_ϵt0, schDa1_qt0.data, schDa1_pS_t0,
    schDa1_ϵt1, pS_tguess )

  return pS_tk1    
end


function update_qn(sch::Schapery, Δt, qnt0, σt0, σt1)
  
  return VectorValue( retqnt1.(
    Ref(sch), sch.λn, 
    Ref(Δt), qnt0.data, 
    Ref(σt0), Ref(σt1)) )
  
end


function stressσ_fnc(seg::Segment, QTr, P, J, ∇u)      
    
	local FΓ, EDir, ETang, stressS, JNew, sΛ
	
	FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
	# FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)
	EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )
	ETang = P ⋅ EDir ⋅ P

	stressS = 2*seg.μm * ETang

	JNew = J  + ∇u'
	sΛ = ((JNew ⊙ JNew) ./ (J ⊙ J)).^0.5

	return ( FΓ ⋅ stressS ⋅ FΓ' ) / sΛ

end


function stressσ_fnc(
  # sch::Schapery, Δt,
  QTr, P, J, ∇u, rotM,
  # schDa1_ϵt0, schDa1_qt0, schDa1_pS_t0, 
  schDa1_pS_t1,
  # schDa2_ϵt0, schDa2_qt0, schDa2_pS_t0, 
  schDa2_pS_t1)
    
	local FΓ, EDir, ETang
  # local pS_tk1, pS_tk2, err1, err2, pS_tguess
  
	
	FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
	# FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)
	# EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )
	# ETang = P ⋅ EDir ⋅ P

  # rotM = getStrRotMatrix( ETang )
  # pETang = rotM ⋅ (ETang ⋅ transpose(rotM))  

  # pS_tguess = schDa1_pS_t0
  # pS_tk1, err1 = get_stressNLVE(
  #   sch, Δt,
  #   schDa1_ϵt0, schDa1_qt0.data, schDa1_pS_t0,
  #   pETang[1], pS_tguess )

  # pS_tguess = schDa2_pS_t0
  # pS_tk2, err2 = get_stressNLVE(
  #   sch, Δt,
  #   schDa2_ϵt0, schDa2_qt0.data, schDa2_pS_t0,
  #   pETang[4], pS_tguess )

  pStr = TensorValue( 
    schDa1_pS_t1,
    0.0, 0.0, 
    schDa2_pS_t1 )

  S = ( transpose(rotM) ⋅ pStr ) ⋅ rotM

	JNew = J  + ∇u'
	sΛ = ((JNew ⊙ JNew) ./ (J ⊙ J)).^0.5

	return ( FΓ ⋅ S ⋅ FΓ' ) / sΛ

end


function ETang_fnc(QTr, P, J, ∇u)      
    
	local FΓ, EDir
	
	FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
	# FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)
	EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )

	return P ⋅ EDir ⋅ P
end
# ----------------------End----------------------



"""
Linear Hooke's Law material
=============

"""
# ---------------------Start---------------------
function linStressStrain(seg, strain::Float64)
  
  ϵErr(σ) = (2*seg.μm*strain - σ)
  σi = find_zero(ϵErr, 2*seg.μm*strain-1e-10)
  
  return σi
end

function linStressStrain(seg::Segment, strain::ForwardDiff.Dual)
  
  # println(typeof(strain))
  return 2*seg.μm * strain
end
# ----------------------End----------------------




"""
Schapery Nonlinear Viscoelastic
=============

Incremental strain formulation

Haj-Ali (2004)

Haj‐Ali, R. M. & Muliana, A. H. 
Numerical finite element formulation of 
the Schapery non‐linear viscoelastic material model. 
Numerical Meth Engineering 59, 25–45 (2004).

"""
# ---------------------Start---------------------


## Aux
retExpTerm1(λn, ΔΨ) = exp(-λn * ΔΨ)

retExpTerm2(λn, ΔΨ) = (1 - exp(-λn * ΔΨ))/(λn * ΔΨ)


## g0 g1 g2
retg0(S::Schapery, σ) = poly_eval(σ, S.g0)
retg1(S::Schapery, σ) = poly_eval(σ, S.g1)
retg2(S::Schapery, σ) = poly_eval(σ, S.g2)


## DBar
function DBar(S::Schapery, ΔΨ, σ)

  local expTerm2

  expTerm2 = retExpTerm2.( S.λn, Ref(ΔΨ) ) 
  # this will be constant if ΔΨ is a constant

  return retg0(S,σ) * S.D0 +
    retg1(S,σ) * retg2(S,σ) * ( 
      S.ΣDn - sum( S.Dn .* expTerm2 ) )

end


## Update qn
function retqnt1(S::Schapery, λn, ΔΨ, qnt0, σt0, σt1)  

  return retExpTerm1(λn, ΔΨ) * qnt0 + 
    ( retg2(S, σt1) * σt1 - retg2(S, σt0) * σt0 ) * 
    retExpTerm2(λn, ΔΨ)
  
end


## Uniaxial stress predictor
function σPredicted( S::Schapery, 
  ϵt0, ΔΨt0, qnt0, σt0,
  ϵt1, ΔΨtk, σtk )
  
  local tmp1, tmp2, tmp3, g1tk, g1t0, g2tk, g2t0
  local DBartk, σtk1, err

  g1t0 = retg1(S, σt0) #known  
  g1tk = retg1(S, σtk) #guess
  g2t0 = retg2(S, σt0) #known  
  g2tk = retg2(S, σtk) #guess
  
  # Aux functions
  tmpfn2( ΔΨt1, g1t0, g1t1, λi, Di, qit0 ) = 
    Di*( g1t1*retExpTerm1(λi, ΔΨt1) - g1t0 )*qit0
  
  tmpfn3(ΔΨt0, ΔΨt1, g1t0, g1t1, λi, Di) = 
    Di*
    ( g1t0*retExpTerm2(λi,ΔΨt0) - g1t1*retExpTerm2(λi,ΔΨt1) )


  tmp1 = DBar(S, ΔΨt0, σt0) * σt0 #known  

  tmp2 = sum( 
    tmpfn2.(Ref(ΔΨtk), Ref(g1t0), Ref(g1tk), S.λn, S.Dn, qnt0) )

  tmp3 = g2t0 * σt0 *
    sum( 
      tmpfn3.(Ref(ΔΨt0), Ref(ΔΨtk), Ref(g1t0), Ref(g1tk), S.λn, S.Dn) 
    )

  DBartk = DBar(S, ΔΨtk, σtk)

  σtk1 = ( ϵt1 - ϵt0 + tmp1 + tmp2 + tmp3 )/DBartk


  ## Calculating residual
  # This is missing the update in ΔΨk for now #TODO

  g1tk = retg1(S, σtk1) #guess
  g2tk = retg2(S, σtk1) #guess

  tmp2 = sum( 
    tmpfn2.(Ref(ΔΨtk), Ref(g1t0), Ref(g1tk), S.λn, S.Dn, qnt0) )

  tmp3 = g2t0 * σt0 *
    sum( 
      tmpfn3.(Ref(ΔΨt0), Ref(ΔΨtk), Ref(g1t0), Ref(g1tk), S.λn, S.Dn) 
    )

  DBartk = DBar(S, ΔΨtk, σtk1)

  err = σtk1 - (ϵt1 - ϵt0 + tmp1 + tmp2 + tmp3)/DBartk
  
  return σtk1, err
end



function Residual_σPredicted( S::Schapery, 
  ϵt0, ΔΨt0, qnt0, σt0,
  ϵt1, ΔΨtk, σtk )
  
  local tmp1, tmp2, tmp3, g1tk, g1t0, g2tk, g2t0
  local DBartk, σtk1, err

  g1t0 = retg1(S, σt0) #known  
  g1tk = retg1(S, σtk) #guess
  g2t0 = retg2(S, σt0) #known  
  g2tk = retg2(S, σtk) #guess
  
  # Aux functions
  tmpfn2( ΔΨt1, g1t0, g1t1, λi, Di, qit0 ) = 
    Di*( g1t1*retExpTerm1(λi, ΔΨt1) - g1t0 )*qit0
  
  tmpfn3(ΔΨt0, ΔΨt1, g1t0, g1t1, λi, Di) = 
    Di*
    ( g1t0*retExpTerm2(λi,ΔΨt0) - g1t1*retExpTerm2(λi,ΔΨt1) )


  tmp1 = DBar(S, ΔΨt0, σt0) * σt0 #known  

  tmp2 = sum( 
    tmpfn2.(Ref(ΔΨtk), Ref(g1t0), Ref(g1tk), S.λn, S.Dn, qnt0) )

  tmp3 = g2t0 * σt0 *
    sum( 
      tmpfn3.(Ref(ΔΨt0), Ref(ΔΨtk), Ref(g1t0), Ref(g1tk), S.λn, S.Dn) 
    )

  DBartk = DBar(S, ΔΨtk, σtk)
  
  err = σtk - ( ϵt1 - ϵt0 + tmp1 + tmp2 + tmp3 )/DBartk

  return abs(err)
end


## Uniaxial viscoelastic strain 
function retϵve( S::Schapery, 
  ϵt0, ΔΨt0, qnt0, σt0,
  ΔΨtk, σtk )
  
  local tmp1, tmp2, tmp3, g1tk, g1t0, g2tk, g2t0
  local DBartk, err

  g1t0 = retg1(S, σt0) #known  
  g1tk = retg1(S, σtk) #guess
  g2t0 = retg2(S, σt0) #known  
  g2tk = retg2(S, σtk) #guess
  
  # Aux functions
  tmpfn2( ΔΨt1, g1t0, g1t1, λi, Di, qit0 ) = 
    Di*( g1t1*retExpTerm1(λi, ΔΨt1) - g1t0 )*qit0
  
  tmpfn3(ΔΨt0, ΔΨt1, g1t0, g1t1, λi, Di) = 
    Di*
    ( g1t0*retExpTerm2(λi,ΔΨt0) - g1t1*retExpTerm2(λi,ΔΨt1) )


  tmp1 = DBar(S, ΔΨt0, σt0) * σt0 #known  

  tmp2 = sum( 
    tmpfn2.(Ref(ΔΨtk), Ref(g1t0), Ref(g1tk), S.λn, S.Dn, qnt0) )

  tmp3 = g2t0 * σt0 *
    sum( 
      tmpfn3.(Ref(ΔΨt0), Ref(ΔΨtk), Ref(g1t0), Ref(g1tk), S.λn, S.Dn) 
    )

  DBartk = DBar(S, ΔΨtk, σtk)

  ϵt1 = ϵt0 - tmp1 - tmp2 - tmp3 + σtk*DBartk


  ## Calculating residual
  # This is missing the update in ΔΨk for now #TODO

  g1tk = retg1(S, σtk) #guess
  g2tk = retg2(S, σtk) #guess

  tmp2 = sum( 
    tmpfn2.(Ref(ΔΨtk), Ref(g1t0), Ref(g1tk), S.λn, S.Dn, qnt0) )

  tmp3 = g2t0 * σt0 *
    sum( 
      tmpfn3.(Ref(ΔΨt0), Ref(ΔΨtk), Ref(g1t0), Ref(g1tk), S.λn, S.Dn) 
    )

  DBartk = DBar(S, ΔΨtk, σtk)

  err = σtk - (ϵt1 - ϵt0 + tmp1 + tmp2 + tmp3)/DBartk
  
  return ϵt1, err
end


# ----------------------End----------------------
end