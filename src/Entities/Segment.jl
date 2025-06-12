module Segments
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.TensorValues
using Gridap.ODEs
using Mooring.Materials: Material
using Mooring.MooringPoints: MooringPoint, MooringPointMotion, get_motion_type

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
  values = Function[]
  for point in s.points
    if get_motion_type(point) !== nothing
      push!(values, get_motion_function(point))
    end
  end
  return values
end

"""
  get_transient_FESpaces(s::Segment, order::Int=1, dim::Int=2)

  Get the transient finite element spaces for a segment. This function creates a test finite element space
  and a transient trial finite element space based on the segment's triangulation, Dirichlet tags, and values.
  It assumes that Lagrangian elements are used.
  Input:
  - `s::Segment`: The segment for which the finite element spaces are created.
  - `order::Int`: The order of the finite element space (default is 1).
  - `dim::Int`: The dimension of the finite element space (default is 2).
Output:
  - `V::TestFESpace`: The test finite element space.
  - `U::TransientTrialFESpace`: The transient trial finite element space.
"""
function get_transient_FESpaces(s::Segment, order::Int=1, dim::Int=2)
  # Get the triangulation of the segment
  trian = get_triangulation(s)

  # Get the Dirichlet tags and values
  dirichlet_tags = get_dirichlet_tags(s)
  dirichlet_values = get_dirichlet_values(s)

  # Create the reference finite element space
  reffe = ReferenceFE(lagrangian, VectorValue{dim, Float64}, order)

  # Create the test finite element space
  V = TestFESpace(trian, reffe, dirichlet_tags=dirichlet_tags)

  # Create the transient trial finite element space
  U = TransientTrialFESpace(V, dirichlet_values)

  return V,U
end

end 