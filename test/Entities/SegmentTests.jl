import Mooring.Segments as Seg
import Mooring.Materials as M
using Gridap

# Create a Discrete model
nelem = 2
model = CartesianDiscreteModel((0,1),(nelem,))
map = (r) -> VectorValue(2r[1],0.0)
material = M.LinearElastic(E=1.0)

# add labels to the model
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"pointA",[1])
add_tag_from_tags!(labels,"pointB",[2])
add_tag_from_tags!(labels,"segment1",[3])

# Create a segment
segment = Seg.Segment(model, "segment1", "pointA", "pointB", map, material)

# Test the segment properties
@test Seg.get_tag(segment) == "segment1"
trian = Seg.get_triangulation(segment)
@test length(get_cell_coordinates(trian)) == nelem
map_func = Seg.get_map(segment)
@test map_func(get_cell_coordinates(trian)[end][end]) == VectorValue(2.0, 0.0)
btrian1, btrian2 = Seg.get_boundary_triangulations(segment)
@test get_cell_coordinates(btrian1)[1][1] == VectorValue(0.0, )
@test get_cell_coordinates(btrian2)[1][1] == VectorValue(1.0,)