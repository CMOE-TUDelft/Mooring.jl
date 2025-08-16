import Mooring.MooringDiscreteModel as DM
import Mooring.Topology as Topo

# Minimal topology
p1 = Topo.TopoPoint(1, [0.0, 0.0], 0.1)
p2 = Topo.TopoPoint(2, [1.0, 0.0], 0.2)
seg = Topo.TopoSegment(1, "line1", p1.id, p2.id, 1.5)

topo = Topo.TopologyData([p1, p2], [seg])

gmsh = DM.generate_mesh(topo)

# Check model name
@test "mooring_model" in gmsh.model.list()

# Retrieve points from gmsh
point_coords = Dict(id => gmsh.model.getValue(0, id, []) for id in [1, 2])
@test isapprox(point_coords[1][1:2], [0.0, 0.0], atol=1e-8)
@test isapprox(point_coords[2][1:2], [1.0, 0.0], atol=1e-8)

# Check that line exists
lines = gmsh.model.getEntities(1) # 1D entities
line_ids = [ent[2] for ent in lines]
@test 1 in line_ids

# Check physical group name
phys_groups = gmsh.model.getPhysicalGroups(1)
@test any(pg -> pg[2] == 1, phys_groups)
name = gmsh.model.getPhysicalName(1, 1)
@test name == "line1"

# Check physical groups for points
phys_groups_points = gmsh.model.getPhysicalGroups(0)
point_ids = [pg[2] for pg in phys_groups_points]
@test all(id in point_ids for id in [1, 2])

# Check physical group names for points
name_p1 = gmsh.model.getPhysicalName(0, 1)
name_p2 = gmsh.model.getPhysicalName(0, 2)
@test name_p1 == "Point_1"
@test name_p2 == "Point_2"

# Check that a mesh was generated
nodes = gmsh.model.mesh.getNodes()
@test length(nodes[1]) > 0  # node ids must exist
@test length(nodes[2]) > 0  # node coords must exist

gmsh.finalize()