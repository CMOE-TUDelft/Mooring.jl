module stressLinear

using Revise
using Gridap
using Parameters



"""
Custom Structs
=============

"""
# ---------------------Start---------------------
struct Seg
  
	ρcDry::Real	
	E::Real
	L::Real
	A::Real

	ρcSub::Real
	μm::Real


  function Seg( ρcDry, E, L, A, ρcSub )

    μm = 0.5*E
    
    new(ρcDry, E, L, A, ρcSub, μm)
  end


	function Seg( params )

		@unpack E, ρcDry, L, A_str, ρw = params
		ρcSub = ρcDry - ρw
    μm = 0.5*E
    
    new(ρcDry, E, L, A_str, ρcSub, μm)
  end
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