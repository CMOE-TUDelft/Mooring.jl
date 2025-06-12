module Segments
using Gridap.Geometry
using Mooring.Materials: Material
using Mooring.MooringPoints: MooringPoint, MooringPointMotion

"""
Segment Struct

This struct is used to define a segment in the mooring system.
It includes the following fields:
- `tag::String`: Identifier for the segment
- `trian::Triangulation`: Triangulation of the segment
- `points::Vector{MooringPoint}`: End points of the segment, defined as [`MooringPoint`](@ref) types.
- `map::Function`: Function to map the segment from reference configuration
to undeformed configuration.
- `material::Material`: [Material](../Physics/Materals) properties of the segment
"""
struct Segment
  tag::String
  trian::Triangulation
  points::Vector{MooringPoint} # boundary triangulations
  map::Function
  material::Material
end

"""
Segment(model::DiscreteModel, segment_tag::String, pointA_tag::String, pointB_tag::String,
        map::Function, material::Material)

Create a segment in the mooring system from a discrete model. The triangulations of the segment 
and boundaries are created from the model. The segment is defined by its tag, the triangulation,
the mapping function, and the material properties.
"""
function Segment(model::DiscreteModel, segment_tag::String, pointA_tag::String, pointB_tag::String,
                 map::Function, material::Material, motionA::MooringPointMotion=nothing, motionB::MooringPointMotion=nothing)

  # Get the triangulation of the segment
  trian = Interior(model, tags=[segment_tag])

  # Define boundary points
  pointA = MooringPoint(model, pointA_tag, motionA)
  pointB = MooringPoint(model, pointB_tag, motionB)

  # Create the segment
  return Segment(segment_tag, trian, [pointA, pointB], map, material)
end

"""
  get_tag(s::Segment)
  Get the tag of a segment
"""
get_tag(s::Segment) = s.tag

"""
  get_triangulation(s::Segment)
  Get the triangulation of a segment
"""
get_triangulation(s::Segment) = s.trian

"""
  get_points(s::Segment)
  Get the boundary points of a segment
"""
get_points(s::Segment) = s.points

"""
  get_map(s::Segment)
  Get the mapping function of a segment. This function maps the segment from the reference configuration
  to the undeformed configuration. In 2D, it is defined as:

  ```math
  X(r) = \begin{bmatrix} f_1(r) \\ f_2(r) \end{bmatrix}
  ```

  where `r` is the reference configuration point.
"""
get_map(s::Segment) = s.map

"""
  get_material(s::Segment)
  Get the material properties of a segment
"""
get_material(s::Segment) = s.material

"""
  get_dirichlet_tags(s::Segment)
  Get the Dirichlet tags of a segment. These tags are used to define the boundary conditions
  for the segment in the finite element space. A segment has as many Dirichlet tags as points with 
  motion different from `nothing`. If both points have motion equal to `nothing`, the segment has no Dirichlet tags.
"""
function get_dirichlet_tags(s::Segment)
  tags = String[]
  for point in s.points
    if get_motion_type(point) !== nothing
      push!(tags, get_tag(point))
    end
  end
  return tags
end

"""
  get_dirichlet_values(s::Segment)
  Get the Dirichlet values of a segment. These values are used to define the boundary conditions
  for the segment in the trial finite element space. A segment has as many Dirichlet values as points with 
  motion different from `nothing`. If both points have motion equal to `nothing`, the segment has no Dirichlet values.
"""
function get_dirichlet_values(s::Segment)
  values = MooringPointMotion[]
  for point in s.points
    if get_motion_type(point) !== nothing
      push!(values, get_motion(point))
    end
  end
  return values
end

end 