module StressNLVE

using Revise
using Gridap
using ForwardDiff
using Parameters
using Roots: find_zero

using Mooring.StressLinear: Segment


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
  g0::Vector{Float64}  
  g1::Vector{Float64}  
  g2::Vector{Float64} 

  # Calculated values
  ΣDn::Real  

end


function Schapery(
  placeholder; #Ask Oriol
  D0::Real,
  Dn::Vector{<:Real} = [0.0], 
  λn::Vector{<:Real} = [1.0],
  g0::Vector{<:Real} = [1.0],
  g1::Vector{<:Real} = [1.0],
  g2::Vector{<:Real} = [1.0] )

  nDn = length(Dn)
  nλn = length(λn)

  N = min(nDn, nλn)

  ΣDn = sum( Dn[1:N] )
  
  Schapery(D0, N, 
    Dn[1:N], λn[1:N], 
    g0, g1, g2, 
    ΣDn)
end



mutable struct SchaperyData
  qt0
  qt1
  ϵt0
  ϵt1
  σt0
  σt1

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
function poly_eval(σ, coeffs::Vector{Float64})
  return sum( c * σ^(i-1) for (i, c) in enumerate(coeffs) )
end
# ----------------------End----------------------



"""
Stress-strain functions
=============

"""
# ---------------------Start---------------------
function stressK_NLVE(seg::Segment, S::Schapery, 
  QTr, P, ∇u)
    
	local FΓ, EDir, ETang
	
	FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
	# FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)
	EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )
	ETang = P ⋅ EDir ⋅ P

	# return 2*seg.μm * (FΓ ⋅ ETang)  

  rotM = getStrRotMatrix( ETang )
  pETang = rotM ⋅ (ETang ⋅ transpose(rotM))  

  pStr = TensorValue( 
    linStressStrain(seg, pETang[1]),
    0.0, 0.0, 
    0.0 )
    # linStressStrain(seg, pETang[4]) )

  S = ( transpose(rotM) ⋅ pStr ) ⋅ rotM

  return FΓ ⋅ S
end


function stressK_damp_fnc(seg::Segment, QTr, P, ∇u, ∇v)
    
	local FΓ, EDirdot, ETangdot, FΓdot
	
	FΓdot = ∇v' ⋅ QTr
	
	FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
	# FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)
	EDirdot = 0.5 * ( FΓdot' ⋅ FΓ + FΓ' ⋅ FΓdot )
	ETangdot = P ⋅ EDirdot ⋅ P

	return 2*seg.μm * seg.c * (FΓ ⋅ ETangdot)
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


## Unixial stress predictor
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

  err = ϵt1 - ϵt0 + tmp1 + tmp2 + tmp3 - DBartk*σtk1
  
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
  
  err = ϵt1 - ϵt0 + tmp1 + tmp2 + tmp3 - DBartk*σtk

  return err
end


# ----------------------End----------------------
end