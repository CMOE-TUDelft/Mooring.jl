module Driver

import Mooring.MooringPoints as Pts
import Mooring.MooringSegments as Seg
import Mooring.Materials as M
import Mooring.TangentialDiffCalculus as TDC
import Mooring.PointMotion as PM
using Gridap
using Gridap.ReferenceFEs: get_order

# Create a Discrete model
nelem = 10
model = CartesianDiscreteModel((0,1),(nelem,))
# add labels to the model
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"pointA",[1])
add_tag_from_tags!(labels,"pointB",[2])
add_tag_from_tags!(labels,"segment1",[3])

# Define the mapping function and material
map = (r) -> VectorValue(2r[1],4r[1]^2)
material = M.LinearElastic(E=1.0e3)

# Create points
pointAmotion_function = (t,x) -> VectorValue(0.0,0.0)
pointBmotion_function = (t,x) -> VectorValue(sin(0.1*t),0.0)
pointAmotion = PM.CustomMotionType(pointAmotion_function)
pointBmotion = PM.CustomMotionType(pointBmotion_function)
pointA = Pts.MooringPoint(model, "pointA", pointAmotion)
pointB = Pts.MooringPoint(model, "pointB", pointBmotion)

# Create a segment
segment = Seg.MooringSegment(model, "segment1", pointA, pointB, map, material)

# Crete FE spaces
V,U = Seg.get_transient_FESpaces(segment)

# # Get Measures
# dΩ, dΓ1, dΓ2 = Seg.get_measures(segment)

# Get reference configuration
Xₕ = Seg.get_reference_configuration(segment, U(0.0))
writevtk(Seg.get_triangulation(segment),"Xh", cellfields=["Xh"=>Xₕ])

# Get quasi-static residual
res = Seg.get_quasi_static_residual(segment, Xₕ)
op = FEOperator(res, U(1.0), V)
nls = NLSolver(LUSolver(), iterations=200, ftol=1e-6, show_trace=true, method=:newton)
uₕ = solve(nls, op)
writevtk(Seg.get_triangulation(segment),"uh", cellfields=["uh"=>uₕ,"uh+Xh"=>uₕ+Xₕ])


  # res0(u, ψu) =          
  #   ∫( ( (∇(ψu)' ⋅ QTrans_cs) ⊙ 
  #     (stressK_fnc∘(QTrans_cs, P_cs, ∇(u))) )*JJ_cs )dΩ +
  #   ∫( ( -ψu ⋅ FWeih_cs )*JJ_cs )dΩ + 
  #   ∫( ( -ψu ⋅ VectorValue(0.0,1.0) * 
  #     (bedSpring_fnc∘(Xh_cs, csTup1...,u, ∇(u), 0.0*u)) )*JJ_cs )dΩ +
  #   ∫( ( 
  #     -ψu ⋅ ( drag_ΓX(0.0, cnstTup1...)
  #       ∘(UCur_cs, waveVel_cs, csTup1...,∇(u), 0.0*u) ) 
  #   )*JJ_cs )dΩ 

end