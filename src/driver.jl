module Driver

import Mooring.Segments as Seg
import Mooring.Materials as M
using Gridap

nelem = 10
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

# Crete FE spaces
order = 1
reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
trian = Seg.get_triangulation(segment)
dirichlet_tags = Seg.get_dirichlet_tags(segment)
dirichlet_values = Seg.get_dirichlet_values(segment)
V = TestFESpace(reffe, trian, dirichlet_tags=dirichlet_tags)


end