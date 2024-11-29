module StressLinear

using Revise
using Gridap
using Parameters

using Mooring.Drag



"""
Custom Structs
=============

"""
# ---------------------Start---------------------
@with_kw struct Segment
  
	ρcDry::Real	
	E::Real
	L::Real
	A::Real	  
  nd::Real # Nominal diameter  
  c::Real # Damping coeff  

  cOnFlag::Bool # Damping on or off ?
	ρcSub::Real
	μm::Real     
  
  dragProp::Drag.DragProperties #Default Drag.NoDrag

end


function Segment( ρcDry, E, L, A, nd, c, ρcSub;
  dragProp = Drag.DragProperties(Drag.NoDrag()) )

  μm = 0.5*E    
  cOnFlag = false
  if(abs(c) > 1e-10)
    cOnFlag = true
  end
  
  Segment(ρcDry, E, L, A, nd, c, 
    cOnFlag, ρcSub, μm, dragProp)
end


function Segment( params )

  @unpack E, ρcDry, L, AStr, ρw, nd, dragProp = params
  @unpack materialDampCoeff = params

  ρcSub = ρcDry - ρw        

  Segment(ρcDry, E, L, AStr, nd, materialDampCoeff, ρcSub;
    dragProp = dragProp )
end
# ----------------------End----------------------



"""
Functions
=============

"""
# ---------------------Start---------------------

# function stressK_fnc(u)
    
# 	local FΓ, EDir, ETang
	
# 	FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)
# 	EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )
# 	ETang = P_cs ⋅ EDir ⋅ P_cs

# 	return 2*μₘ * (FΓ ⋅ ETang)
# end


function stressK_fnc(seg, QTr, P, ∇u)
    
	local FΓ, EDir, ETang
	
	FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
	# FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)
	EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )
	ETang = P ⋅ EDir ⋅ P

	return 2*seg.μm * (FΓ ⋅ ETang)
end


# function stressK_damp_fnc(seg, QTr, P, v)
    
# 	local FΓ, EDir, ETang, ∇v, tmpTens

#   ∇v = v / seg.L # approximation of strain rate

#   tmpTens = TensorValue(
#     ∇v[1]*QTr[1], ∇v[1]*QTr[2], 
#     ∇v[2]*QTr[1], ∇v[2]*QTr[2]
#   )
	
# 	# FΓ = ( ∇v' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
#   FΓ = tmpTens + TensorValue(1.0,0.0,0.0,1.0)
# 	EDir = 0.5 * ( FΓ' ⋅ FΓ - TensorValue(1.0,0.0,0.0,1.0) )
# 	ETang = P ⋅ EDir ⋅ P

# 	return 2*seg.μm * seg.c * (FΓ ⋅ ETang)
# end

function stressK_damp_fnc(seg, QTr, P, ∇u, ∇v)
    
	local FΓ, EDirdot, ETangdot, FΓdot
	
	FΓdot = ∇v' ⋅ QTr
	
	FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
	# FΓ = ∇(u)' ⋅ QTrans_cs + TensorValue(1.0,0.0,0.0,1.0)
	EDirdot = 0.5 * ( FΓdot' ⋅ FΓ + FΓ' ⋅ FΓdot )
	ETangdot = P ⋅ EDirdot ⋅ P

	return 2*seg.μm * seg.c * (FΓ ⋅ ETangdot)
end


function stressσ_fnc(seg, QTr, P, J, ∇u)      
    
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
    
end