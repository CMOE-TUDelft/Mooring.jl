import Mooring.Segments as Seg
import Mooring.MooringPoints as Pt
import Mooring.Materials as M
using Gridap
using Gridap.FESpaces: SingleFieldFESpace, TransientTrialFESpace

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
point1, point2 = Seg.get_points(segment)
btrian1 = Pt.get_triangulation(point1)
btrian2 = Pt.get_triangulation(point2)
@test get_cell_coordinates(btrian1)[1][1] == VectorValue(0.0, )
@test get_cell_coordinates(btrian2)[1][1] == VectorValue(1.0,)

# Test segment FE spaces
dirichlet_tags = Seg.get_dirichlet_tags(segment)
@test dirichlet_tags == String[]
dirichlet_values = Seg.get_dirichlet_values(segment)
@test dirichlet_values == Function[]
V, U = Seg.get_transient_FESpaces(segment)
@test isa(V, SingleFieldFESpace)
@test isa(U, TransientTrialFESpace)

# Test segment measures
dΩ, dΓ1, dΓ2 = Seg.get_measures(segment)
@test isa(dΩ, Measure)
@test isa(dΓ1, Measure)
@test isa(dΓ2, Measure)