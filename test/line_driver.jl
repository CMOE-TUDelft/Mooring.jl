module LineDriver

import Mooring.MooringPoints as Pts
import Mooring.MooringSegments as Seg
import Mooring.Materials as M
import Mooring.TangentialDiffCalculus as TDC
import Mooring.PointMotion as PM
import Mooring.MooringTopology as Topo
import Mooring.MooringDiscreteModel as DM
using Gridap
using Gridap.ReferenceFEs: get_order
using GridapGmsh: gmsh, GmshDiscreteModel

# MooringTopology
p1 = Topo.TopoPoint(1, [0.0, 0.0], 1.0)
p2 = Topo.TopoPoint(2, [5.0, 3.0], 1.0)
p3 = Topo.TopoPoint(3, [10.0, 5.0], 1.0)
s1 = Topo.TopoSegment(1, 1, 2, 10.0)
s2 = Topo.TopoSegment(2, 2, 3, 15.0)
topo = Topo.MooringTopologyData([p1,p2,p3],[s1,s2])

# Discrete model
model = DM.generate_discrete_model(topo)
writevtk(model, "model_geo_points_lines")

# Define the mapping function and material
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

# Get quasi-static residual
res1 = Seg.get_quasi_static_residual(segment1, X1ₕ)
op1 = FEOperator(res1, U1(1.0), V1)
nls = NLSolver(LUSolver(), iterations=200, ftol=1e-6, show_trace=true, method=:newton)
u1ₕ = solve(nls, op1)
writevtk(Seg.get_triangulation(segment1),"u1h", cellfields=["uh"=>u1ₕ,"uh+Xh"=>u1ₕ+X1ₕ])

end