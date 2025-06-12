module MooringPoints
using Gridap.Geometry
using Mooring.PointMotion: MotionType, get_point_motion_function

MooringPointMotion = Union{MotionType, Nothing}

"""
MooringPoint Struct
This struct is used to define a point in the mooring system.
It includes the following fields:
- `tag::String`: Identifier for the point
- `btrian::Triangulation`: BoundaryTriangulation or InterfaceTriangulation with the position of the point in space
- `motion::MotionType`: Type of motion of the point. If no motion is specified, the point is free to move in space.
"""
struct MooringPoint
  tag::String
  btrian::Triangulation 
  motion::MooringPointMotion
end

"""
MooringPoint(model::DiscreteModel, point_tag::String, motion::MotionType)
Create a point in the mooring system from a discrete model. The triangulation of the point
and boundaries are created from the model. The point is defined by its tag, the triangulation,
and the motion type (Nothing by default).
"""
function MooringPoint(model::DiscreteModel, point_tag::String, motion::MooringPointMotion=nothing)
  # Get the triangulation of the point
  btrian = Boundary(model, tags=[point_tag])
  
  # Create the point
  return MooringPoint(point_tag, btrian, motion)
end

"""
  get_tag(p::MooringPoint)
  Get the tag of a point
"""
get_tag(p::MooringPoint) = p.tag

"""
  get_triangulation(p::MooringPoint)
  Get the boundary/interface triangulation of a point
"""
get_triangulation(p::MooringPoint) = p.btrian

"""
  get_motion_type(p::MooringPoint)
  Get the motion type of a point
"""
get_motion_type(p::MooringPoint) = p.motion

"""
  get_motion(p::MooringPoint)
  Get the motion function of a point
"""
get_motion_function(p::MooringPoint) = get_point_motion_function(p.motion)

end