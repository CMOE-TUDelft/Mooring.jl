import Mooring.MooringDiscreteModel as DM
import Mooring.ParameterHandler as PH
using GridapGmsh: GmshDiscreteModel
using Gridap.Geometry
using LinearAlgebra

@testset "Topology tests" begin
  ph = PH.ParameterHandler()

  # --- Create sample topology ---
  # 1 --(10.0)-- 2 --(5.0)-- 3
  #              |
  #            (7.0)
  #              |
  #              4
  ph.points[1] = PH.PointParameters(id=1, coords=[1.0, 2.0, 3.0])
  ph.points[2] = PH.PointParameters(id=2, coords=[4.0, 6.0, 8.0])
  ph.points[3] = PH.PointParameters(id=3, coords=[0.0, 0.0, 0.0])
  ph.points[4] = PH.PointParameters(id=4, coords=[0.0, 0.0, 0.0])

  ph.segments[1] = PH.SegmentParameters(id=1, start_point=1, stop_point=2, length=10.0)
  ph.segments[2] = PH.SegmentParameters(id=2, start_point=2, stop_point=3, length=5.0)
  ph.segments[3] = PH.SegmentParameters(id=3, start_point=2, stop_point=4, length=7.0)

  ph.lines[1] = PH.LineParameters(id=1, points=[1,2,3,4], segments=[1,2,3])

  # --- Test build_adjacency ---
  adj = DM.build_adjacency(ph.lines[1],ph)
  @test length(adj) == 4
  @test adj[1] == [(2, 10.0)]
  @test adj[2] == [(1, 10.0), (3, 5.0), (4, 7.0)]
  @test adj[3] == [(2, 5.0)]
  @test adj[4] == [(2, 7.0)]

  # --- Test assign_coords ---
  coords = DM.assign_coords(ph.lines[1], ph; anchor=1)
  @test coords[1] == 0.0
  @test coords[2] == 10.0
  @test coords[3] == 15.0
  @test coords[4] == 17.0

  # Check different anchor
  coords2 = DM.assign_coords(ph.lines[1], ph; anchor=3)
  @test coords2[3] == 0.0
  @test coords2[2] == 5.0
  @test coords2[1] == 15.0
  @test coords2[4] == 12.0

  # # Check map
  # x_p1 = ph.points[1].coords
  # x_p2 = ph.points[2].coords
  # length1 = norm(x_p2 .- x_p1)
  # map1 = Topo.get_physical_map(seg1,topo)

  # @test map1(VectorValue(0.0)) ≈ VectorValue(x_p1)
  # @test map1(VectorValue(length1)) ≈ VectorValue(x_p2)
  # mid1 = map1(VectorValue(length1/2))
  # @test mid1 ≈ VectorValue([(x_p1[i] + x_p2[i])/2 for i in 1:3])

end # end of testset

@testset "Discrete model tests" begin
  ph = PH.ParameterHandler()

  # Minimal topology
  ph.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0], mesh_size=0.1)
  ph.points[2] = PH.PointParameters(id=2, coords=[1.0, 0.0], mesh_size=0.2)
  ph.segments[1] = PH.SegmentParameters(id=1, tag="line1", start_point=1, stop_point=2, length=1.5)
  ph.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])

  gmsh = DM.generate_mesh(ph.lines[1],ph)

  # Check model name
  @test "mooring_model" in gmsh.model.list()

  # Retrieve points from gmsh
  point_coords = Dict(id => gmsh.model.getValue(0, id, []) for id in [1, 2])
  @test isapprox(point_coords[1][1:2], [0.0, 0.0], atol=1e-8)
  @test isapprox(point_coords[2][1:2], [1.5, 0.0], atol=1e-8)

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

  # Test discrete model
  model = DM.generate_discrete_model(ph.lines[1],ph)

  # 1. Check type
  @test model isa UnstructuredDiscreteModel

  # 2. Check topology (entities in model)
  labels = get_face_labeling(model)
  dim_ids = labels.d_to_dface_to_entity
  point_ids = dim_ids[1]
  line_ids = dim_ids[2]
  id_to_name = labels.tag_to_name
  @test id_to_name[point_ids[1]] == "Point_1"
  @test id_to_name[point_ids[2]] == "Point_2"
  @test id_to_name[line_ids[1]] == "line1"

end # end of testset