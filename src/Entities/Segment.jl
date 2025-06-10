module Segments
using Gridap.Geometry
using Mooring.Materials: Material

"""
Segment Struct

This struct is used to define a segment in the mooring system.
It includes the following fields:
- `tag::String`: Identifier for the segment
- `trian::Triangulation`: Triangulation of the segment
- `btrians::Vector{Triangulation}`: Boundary triangulations of the segment
- `map::Function`: Function to map the segment from reference configuration
to undeformed configuration.
- `material::Material`: [Material](../Physics/Materals) properties of the segment
"""
struct Segment
  tag::String
  trian::Triangulation
  btrians::Vector{Triangulation} # boundary triangulations
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
                 map::Function, material::Material)

  # Get the triangulation of the segment
  trian = Interior(model, tags=[segment_tag])
  btrians = [Boundary(model, tags=[pointA_tag]), Boundary(model, tags=[pointB_tag])]
  
  # Create the segment
  return Segment(segment_tag, trian, btrians, map, material)
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
  get_boundary_triangulations(s::Segment)
  Get the boundary triangulations of a segment
"""
get_boundary_triangulations(s::Segment) = s.btrians

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


end 