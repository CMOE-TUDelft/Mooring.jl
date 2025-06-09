module Segment

"""
Segment Struct

This struct is used to define a segment in the mooring system.
It includes the following fields:
- `tag::String`: Identifier for the segment
- `trian::Triangulation`: Triangulation of the segment
- `map::Function`: Function to map the segment from reference configuration
to undeformed configuration.
- `material::Material`: [Material](../Physics/Materals) properties of the segment
"""
struct Segment
  tag::String
  trian::Triangulation
  map::Function
  material::Material
end

end 