import Mooring.Topology as Topo

# --- Create sample topology ---
# 1 --(10.0)-- 2 --(5.0)-- 3
#              |
#            (7.0)
#              |
#              4
p1 = Topo.TopoPoint(1, [0.0, 0.0, 0.0])
p2 = Topo.TopoPoint(2, [0.0, 0.0, 0.0])
p3 = Topo.TopoPoint(3, [0.0, 0.0, 0.0])
p4 = Topo.TopoPoint(4, [0.0, 0.0, 0.0])

seg1 = Topo.TopoSegment(1, 1, 2, 10.0, r -> (r, 0.0, 0.0))
seg2 = Topo.TopoSegment(2, 2, 3, 5.0,  r -> (r, 1.0, 0.0))
seg3 = Topo.TopoSegment(3, 2, 4, 7.0,  r -> (r, 3.0, 0.0))

topo = Topo.TopologyData([p1, p2, p3, p4], [seg1, seg2, seg3])

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