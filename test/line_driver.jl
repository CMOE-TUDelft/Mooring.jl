module LineDriver

import Mooring.MooringPoints as Pts
import Mooring.MooringSegments as Seg
import Mooring.Materials as M
import Mooring.TangentialDiffCalculus as TDC
import Mooring.PointMotion as PM
import Mooring.Topology as Topo
import Mooring.MooringDiscreteModel as DM
using Gridap
using Gridap.ReferenceFEs: get_order
using GridapGmsh: gmsh, GmshDiscreteModel

# # Create Gmsh Geometry
# gmsh.initialize()
# lc = 1e-1  # mesh size
# p1x = 0.0; p1y = 0.0; p1z = 0.0
# p2x = 0.0; p2y = 1.0; p2z = 0.0
# p3x = 10.0; p3y = 5.0; p3z = 0.0

# # Define 3 points
# p1 = gmsh.model.geo.addPoint(p1x, p1y, p1z, lc, 1)
# p2 = gmsh.model.geo.addPoint(p2x, p2y, p2z, lc, 2)
# p3 = gmsh.model.geo.addPoint(p3x, p3y, p3z, lc, 3)

# # Define 2 straight lines: line1 from p1→p2, line2 from p2→p3
# l1 = gmsh.model.geo.addLine(p1, p2, 1)
# l2 = gmsh.model.geo.addLine(p2, p3, 2)

# # Synchronize geometry
# gmsh.model.geo.synchronize()

# # Add physical groups for each point (dim=0)
# pg_p1 = gmsh.model.addPhysicalGroup(0, [p1], 11)
# gmsh.model.setPhysicalName(0, pg_p1, "PointA")
# pg_p2 = gmsh.model.addPhysicalGroup(0, [p2], 12)
# gmsh.model.setPhysicalName(0, pg_p2, "PointB")
# pg_p3 = gmsh.model.addPhysicalGroup(0, [p3], 13)
# gmsh.model.setPhysicalName(0, pg_p3, "PointC")

# # Add physical groups for each line (dim=1)
# pg_l1 = gmsh.model.addPhysicalGroup(1, [l1], 21)
# gmsh.model.setPhysicalName(1, pg_l1, "Line1")
# pg_l2 = gmsh.model.addPhysicalGroup(1, [l2], 22)
# gmsh.model.setPhysicalName(1, pg_l2, "Line2")

# # (Optional) Generate mesh
# gmsh.model.mesh.generate(1)

# # Optional: Write mesh to file
# gmsh.write("model_geo_points_lines.msh")

# # Finalize Gmsh session
# gmsh.finalize()

# # Create a Discrete model
# model = GmshDiscreteModel("model_geo_points_lines.msh")
# writevtk(model, "model_geo_points_lines.vtu")

# Topology
p1 = Topo.TopoPoint(1, [0.0, 0.0], 1.0)
p2 = Topo.TopoPoint(2, [5.0, 3.0], 1.0)
p3 = Topo.TopoPoint(3, [10.0, 5.0], 1.0)
s1 = Topo.TopoSegment(1, 1, 2, 10.0)
s2 = Topo.TopoSegment(2, 2, 3, 15.0)
topo = Topo.TopologyData([p1,p2,p3],[s1,s2])

# Discrete model
model = DM.generate_discrete_model(topo)
writevtk(model, "model_geo_points_lines")

# Define the mapping function and material
# function map(r)
#   fx = p1x + (p2x - p1x) * (r[1]-p1x) / (p2x - p1x)
#   yc = (p1x+p2x)/2
#   h = (p2y - p1y) / 2
#   fy = 1+((r[1]-yc)/h)
#   println("r: ", r, "map: ",VectorValue(r[1], fx*fy))  
#   VectorValue(r[1], fx*fy)  
# end
map1 = Topo.get_physical_map(s1, topo)
map2 = Topo.get_physical_map(s2, topo)
material1 = M.LinearElastic(E=1.0e3)
material2 = M.LinearElastic(E=1.0e0)

# Create points
p1x,p1y = Topo.get_coords(p1)
p2x,p2y = Topo.get_coords(p2)
p3x,p3y = Topo.get_coords(p3)
pointAmotion_function = (t,x) -> VectorValue(0.0,0.0)#ß(p1x,p1y)
pointBmotion_function = (t,x) -> VectorValue(p2x+sin(0.1*t),p2y)
pointCmotion_function = (t,x) -> VectorValue(p3x+sin(0.5*t),p3y)
pointAmotion = PM.CustomMotionType(pointAmotion_function)
pointBmotion = PM.CustomMotionType(pointBmotion_function)
pointCmotion = PM.CustomMotionType(pointCmotion_function)
pointA = Pts.MooringPoint(model, p1.tag, pointAmotion)
pointB = Pts.MooringPoint(model, p2.tag, pointBmotion)
pointC = Pts.MooringPoint(model, p3.tag, pointCmotion)

# Create a segments
segment1 = Seg.MooringSegment(model, s1.tag, pointA, pointB, map1, material1, 10.0)
segment2 = Seg.MooringSegment(model, s2.tag, pointB, pointC, map2, material2, 1.0)

# Create FE spaces
V1,U1 = Seg.get_transient_FESpaces(segment1)
V2,U2 = Seg.get_transient_FESpaces(segment2)

# # Get Measures
dΩ1, dΓ1, dΓ2 = Seg.get_measures(segment1)

# Get reference configuration
X1ₕ = Seg.get_reference_configuration(segment1, U1(0.0))
X2ₕ = Seg.get_reference_configuration(segment2, U2(0.0))
writevtk(Seg.get_triangulation(segment1), "X1h", cellfields=["X1h"=>X1ₕ])
writevtk(Seg.get_triangulation(segment2), "X2h", cellfields=["X2h"=>X2ₕ])

# TDC quantities
J = TDC.J(X1ₕ)
Jabs = norm∘(J)
Q = TDC.Q∘J

# Xq = get_cell_points(dΩ1)
# println(Jabs(Xq))

# Define the stress functions
S(u) = M.S(material1, X1ₕ, u)
FΓ(u) = TDC.FΓ(u, X1ₕ)
K(u) = M.K∘(FΓ(u), S(u))

# Define gravity force
Fᵨ = VectorValue(0.0, -Seg.get_density(segment1) * Seg.get_area(segment1) * 9.81)


# Get quasi-static residual
res1 = Seg.get_quasi_static_residual(segment1, X1ₕ)
# res1(u,v) = ∫((∇(u)⊙∇(v))*Jabs)dΩ1
# res1(u,v) = ∫(((∇(v)' ) ⊙ S(u)) * Jabs)dΩ1
op1 = FEOperator(res1, U1(1.0), V1)
nls = NLSolver(LUSolver(), iterations=200, ftol=1e-6, show_trace=true, method=:newton)
u1ₕ = solve(nls, op1)
writevtk(Seg.get_triangulation(segment1),"u1h", cellfields=["uh"=>u1ₕ,"uh+Xh"=>u1ₕ+X1ₕ])

#=
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

  =#

end