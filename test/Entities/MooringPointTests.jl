import Mooring.MooringPoints as Pt
using Gridap

# Create a Discrete model
nelem = 2
model = CartesianDiscreteModel((0,1),(nelem,))

# add labels to the model
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"pointA",[1])
add_tag_from_tags!(labels,"pointB",[2])

# Create triangulations
trian = Triangulation(model)

# Create points
s_id_trians = [(1, trian)]
pointA = Pt.MooringPoint(s_id_trians, "pointA", nothing)
pointB = Pt.MooringPoint(s_id_trians, "pointB")

# Test the point properties
@test Pt.get_tag(pointA) == "pointA"
@test Pt.get_tag(pointB) == "pointB"
trianA = Pt.get_triangulation(pointA)
trianB = Pt.get_triangulation(pointB)
@test get_cell_coordinates(trianA)[1][1] == VectorValue(0.0, )
@test get_cell_coordinates(trianB)[1][1] == VectorValue(1.0,)
motionA = Pt.get_motion_function(pointA)
@test motionA == nothing
motionB = Pt.get_motion_function(pointB)
@test motionB == nothing