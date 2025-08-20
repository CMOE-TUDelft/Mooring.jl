import Mooring.MooringTopology as Topo
using LinearAlgebra

# --- Create sample topology ---
# 1 --(10.0)-- 2 --(5.0)-- 3
#              |
#            (7.0)
#              |
#              4
p1 = Topo.TopoPoint(1, [1.0, 2.0, 3.0])
p2 = Topo.TopoPoint(2, [4.0, 6.0, 8.0])
p3 = Topo.TopoPoint(3, [0.0, 0.0, 0.0])
p4 = Topo.TopoPoint(4, [0.0, 0.0, 0.0])

seg1 = Topo.TopoSegment(1, 1, 2, 10.0)
seg2 = Topo.TopoSegment(2, 2, 3, 5.0)
seg3 = Topo.TopoSegment(3, 2, 4, 7.0)

topo = Topo.MooringTopologyData([p1, p2, p3, p4], [seg1, seg2, seg3])

# --- Test build_adjacency ---
adj = Topo.build_adjacency(topo)
@test length(adj) == 4
@test adj[1] == [(2, 10.0)]
@test adj[2] == [(1, 10.0), (3, 5.0), (4, 7.0)]
@test adj[3] == [(2, 5.0)]
@test adj[4] == [(2, 7.0)]

# --- Test assign_coords ---
coords = Topo.assign_coords(topo; anchor=1)
@test coords[1] == 0.0
@test coords[2] == 10.0
@test coords[3] == 15.0
@test coords[4] == 17.0

# Check different anchor
coords2 = Topo.assign_coords(topo; anchor=3)
@test coords2[3] == 0.0
@test coords2[2] == 5.0
@test coords2[1] == 15.0
@test coords2[4] == 12.0

# Check map
x_p1 = Topo.get_coords(p1)
x_p2 = Topo.get_coords(p2)
length1 = norm(x_p2 .- x_p1)
map1 = Topo.get_physical_map(seg1,topo)

@test map1(VectorValue(0.0)) ≈ VectorValue(x_p1)
@test map1(VectorValue(length1)) ≈ VectorValue(x_p2)
mid1 = map1(VectorValue(length1/2))
@test mid1 ≈ VectorValue([(x_p1[i] + x_p2[i])/2 for i in 1:3])