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
  nd::Real

	ρcSub::Real
	μm::Real   
  
  dragProp::Drag.DragProperties #Default Drag.NoDrag

end


function Segment( ρcDry, E, L, A, nd, ρcSub;
  dragProp = Drag.DragProperties(Drag.NoDrag()) )

  μm = 0.5*E    
  
  Segment(ρcDry, E, L, A, nd, ρcSub, μm, dragProp)
end


function Segment( params )

  @unpack E, ρcDry, L, A_str, ρw, nd, dragType = params

  ρcSub = ρcDry - ρw    

  if( typeof(dragType) == Drag.NoDrag )    
    return Segment(ρcDry, E, L, A_str, nd, ρcSub )
  end  

  dragProp = Drag.DragProperties(
    dragType, ρw, nd, A_str
  )

  Segment(ρcDry, E, L, A_str, nd, ρcSub;
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