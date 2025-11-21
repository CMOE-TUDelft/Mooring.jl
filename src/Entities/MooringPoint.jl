module MooringPoints
using Gridap.Geometry
using Gridap.Geometry: get_faces, get_node_coordinates, CompositeTriangulation
using Mooring.PointMotion: MotionType, get_point_motion_function

MooringPointMotion = Union{MotionType, Nothing}

"""
MooringPoint Struct
This struct is used to define a point in the mooring system.
It includes the following fields:
- `tag::String`: Identifier for the point
- `btrian::Triangulation`: BoundaryTriangulation or InterfaceTriangulation with the position of the point in space
- `strian::Vector{Tuple{Int,Triangulation}}`: Vector of tuples containing segment ID and corresponding Triangulation connected to this point
- `motion::MotionType`: Type of motion of the point. If no motion is specified, the point is free to move in space.
"""
struct MooringPoint
  tag::String
  btrian::Triangulation 
  s_id_trians::Vector{Tuple{Int,<:Triangulation}}
  motion::MooringPointMotion
end

"""
MooringPoint(model::DiscreteModel, point_tag::String, motion::MotionType)
Create a point in the mooring system from a discrete model. The triangulation of the point
and boundaries are created from the model. The point is defined by its tag, the triangulation,
and the motion type (Nothing by default).
"""
function MooringPoint(s_id_trians, point_tag::String, motion::MooringPointMotion=nothing)

  # Get the triangulation of the point
  if length(s_id_trians) == 1
    trian = s_id_trians[1][2]
    btrian = Boundary(trian, tags=[point_tag])
  elseif length(s_id_trians) == 2
    trian1 = s_id_trians[1][2]
    trian2 = s_id_trians[2][2]
    btrian = Interface(trian1, trian2)
  else
    error("No cells attached to point with tag $point_tag, or the point is attached to more than two cells (not supported yet).")
  end
 
  # Create the point
  return MooringPoint(point_tag, btrian, s_id_trians, motion)
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

"""
  get_background_node_from_trian(trian::CompositeTriangulation/SkeletonTriangulation)
  Get the background node ID from a point part of a CompositeTriangulation (resulting from a boundary triangulation)
  or a SkeletonTriangulation (resulting from an interface triangulation)
"""
get_background_node_from_trian(trian::CompositeTriangulation) = trian.dtrian.glue.face_to_bgface[1]
get_background_node_from_trian(trian::SkeletonTriangulation) = trian.plus.glue.face_to_bgface[1]

"""
  get_reference_node_coord(p::MooringPoint)
  Get the node coordinates of a point in the background triangulation
"""
function get_reference_node_coord(p::MooringPoint)
  p_bg_id = get_background_node_from_trian(p.btrian)
  coords = get_node_coordinates(p.btrian)[p_bg_id][1]
  return coords
end

end