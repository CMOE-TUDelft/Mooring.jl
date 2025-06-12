module Driver

# import Mooring.MooringPoints as Pts
import Mooring.Segments as Seg
import Mooring.Materials as M
import Mooring.TangentialDiffCalculus as TDC
import Mooring.PointMotion as PM
using Gridap
using Gridap.ReferenceFEs: get_order

nelem = 10
model = CartesianDiscreteModel((0,1),(nelem,))
map = (r) -> VectorValue(2r[1],4r[1]^2)
material = M.LinearElastic(E=1.0)

# add labels to the model
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"pointA",[1])
add_tag_from_tags!(labels,"pointB",[2])
add_tag_from_tags!(labels,"segment1",[3])

# Create a segment
pointBmotion_function = (t,x) -> VectorValue(sin(0.1*t),0.0)
segment = Seg.Segment(model, "segment1", "pointA", "pointB", map, material, nothing, PM.CustomMotionType(pointBmotion_function))

# Crete FE spaces
V,U = Seg.get_transient_FESpaces(segment)

# Get Measures
dΩ, dΓ1, dΓ2 = Seg.get_measures(segment)

# Get reference configuration
Xₕ = Seg.get_reference_configuration(segment, U(0.0))
J = TDC.J(Xₕ)
Jabs = norm∘(J)
Q = TDC.Q∘J

material = Seg.get_material(segment)
S(u) = M.S(material,Xₕ,u)
FΓ(u) = TDC.FΓ(u,Xₕ)
K(u) = M.K∘(FΓ(u),S(u))
writevtk(Seg.get_triangulation(segment),"Xh", cellfields=["Xh"=>Xₕ,"Q"=>Q, "J"=>J, "Jabs"=>Jabs])
res(u,v) = ∫( ((∇(v)'⋅Q') ⊙ K(u) ) * Jabs )dΩ
op = FEOperator(res, U(1.0), V)
uₕ = solve(op)
writevtk(Seg.get_triangulation(segment),"uh", cellfields=["uh"=>uₕ])


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