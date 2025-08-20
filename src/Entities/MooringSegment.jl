module MooringSegments
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.TensorValues
using Gridap.ODEs
using Gridap.CellData
import Mooring.Materials as M
import Mooring.MooringPoints as Pts
import Mooring.TangentialDiffCalculus as TDC

"""
MooringSegment Struct

This struct is used to define a segment in the mooring system.
It includes the following fields:
- `tag::String`: Identifier for the segment
- `trian::Triangulation`: Triangulation of the segment
- `points::Vector{MooringPoint}`: End points of the segment, defined as [`MooringPoint`](@ref) types.
- `map::Function`: Function to map the segment from reference configuration
to undeformed configuration.
- `material::Material`: [Material](../Physics/Materals) properties of the segment
- `density::Real`: Submerged density of the segment (default is 1.0)
- `area::Real`: Effective cross-sectional area of the segment (default is 1.0)
"""
struct MooringSegment
  tag::String
  trian::Triangulation
  points::Vector{Pts.MooringPoint} # boundary triangulations
  map::Function
  material::M.Material
  density::Real
  area::Real
end

"""
MooringSegment(model::DiscreteModel, segment_tag::String, pointA::Pts.MooringPoint, pointB::Pts.MooringPoint,
        map::Function, material::Material, density::Real=1.0, area::Real=1.0)

Create a segment in the mooring system from a discrete model. The triangulations of the segment 
and boundaries are created from the model. The segment is defined by its tag, the triangulation,
the mapping function, and the material properties.
"""
function MooringSegment(model::DiscreteModel, 
                 segment_tag::String, 
                 pointA::Pts.MooringPoint, 
                 pointB::Pts.MooringPoint,
                 map::Function, 
                 material::M.Material,
                 density::Real=1.0, 
                 area::Real=1.0)  

  # Get the triangulation of the segment
  trian = Interior(model, tags=[segment_tag])

  # Create the segment
  return MooringSegment(segment_tag, trian, [pointA, pointB], map, material, density, area)
end

"""
  get_tag(s::MooringSegment)
  Get the tag of a segment
"""
get_tag(s::MooringSegment) = s.tag

"""
  get_triangulation(s::MooringSegment)
  Get the triangulation of a segment
"""
get_triangulation(s::MooringSegment) = s.trian

"""
  get_points(s::MooringSegment)
  Get the boundary points of a segment
"""
get_points(s::MooringSegment) = s.points

"""
  get_map(s::MooringSegment)
  Get the mapping function of a segment. This function maps the segment from the reference configuration
  to the undeformed configuration. In 2D, it is defined as:

  ```math
  X(r) = \begin{bmatrix} f_1(r) \\ f_2(r) \end{bmatrix}
  ```

  where `r` is the reference configuration point.
"""
get_map(s::MooringSegment) = s.map

"""
  get_material(s::MooringSegment)
  Get the material properties of a segment
"""
get_material(s::MooringSegment) = s.material

"""
  get_density(s::MooringSegment)
  Get the density of a segment
"""
get_density(s::MooringSegment) = s.density

"""
  get_area(s::MooringSegment)
  Get the area of a segment
"""
get_area(s::MooringSegment) = s.area

"""
  get_dirichlet_tags(s::MooringSegment)
  Get the Dirichlet tags of a segment. These tags are used to define the boundary conditions
  for the segment in the finite element space. A segment has as many Dirichlet tags as points with 
  motion different from `nothing`. If both points have motion equal to `nothing`, the segment has no Dirichlet tags.
"""
function get_dirichlet_tags(s::MooringSegment)
  tags = String[]
  for point in s.points
    if Pts.get_motion_type(point) !== nothing
      push!(tags, Pts.get_tag(point))
    end
  end
  return tags
end

"""
  get_dirichlet_values(s::MooringSegment)
  Get the Dirichlet values of a segment. These values are used to define the boundary conditions
  for the segment in the trial finite element space. A segment has as many Dirichlet values as points with 
  motion different from `nothing`. If both points have motion equal to `nothing`, the segment has no Dirichlet values.
"""
function get_dirichlet_values(s::MooringSegment)
  values = Function[]
  for point in s.points
    if Pts.get_motion_type(point) !== nothing
      push!(values, Pts.get_motion_function(point))
    end
  end
  return values
end

"""
  get_transient_FESpaces(s::MooringSegment, order::Int=1, dim::Int=2)

  Get the transient finite element spaces for a segment. This function creates a test finite element space
  and a transient trial finite element space based on the segment's triangulation, Dirichlet tags, and values.
  It assumes that Lagrangian elements are used.
  Input:
  - `s::MooringSegment`: The segment for which the finite element spaces are created.
  - `order::Int`: The order of the finite element space (default is 1).
  - `dim::Int`: The dimension of the finite element space (default is 2).
Output:
  - `V::TestFESpace`: The test finite element space.
  - `U::TransientTrialFESpace`: The transient trial finite element space.
"""
function get_transient_FESpaces(s::MooringSegment, order::Int=1, dim::Int=2)
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

"""
  get_measures(s::MooringSegment, degree::Int=2)

  Get the measures for a segment. This function returns the measure for the segment and the boundary points.
  Input:
  - `s::MooringSegment`: The segment for which the measures are created.
  - `degree::Int`: The degree of the measure (default is 2).
Output:
  - `dΩ::Measure`: The measure for the segment.
  - `dΓ1::Measure`: The measure for the first boundary point.
  - `dΓ2::Measure`: The measure for the second boundary point.
"""
function get_measures(s::MooringSegment, degree::Int=2)
  # Get the triangulation of the segment
  trian = get_triangulation(s)

  # Create the measure for the segment
  dΩ = Measure(trian, degree)

  # Get the boundary points
  point1, point2 = get_points(s)

  # Create measures for the boundary points
  dΓ1 = Measure(Pts.get_triangulation(point1), degree)
  dΓ2 = Measure(Pts.get_triangulation(point2), degree)

  return dΩ, dΓ1, dΓ2
end

"""
  get_reference_configuration(s::MooringSegment, U::SingleFieldFESpace)

  Get the reference configuration of a segment. This function uses the map function of the segment
  to create a FEFunction that represents the reference configuration.
  Input:
  - `s::MooringSegment`: The segment for which the reference configuration is created.
  - `U::SingleFieldFESpace`: The finite element space used for interpolation.
Output:
  - `Xₕ::FEFunction`: The FEFunction representing the reference configuration of the segment.
"""
function get_reference_configuration(s::MooringSegment, U::SingleFieldFESpace)
  # Get the map function of the segment
  map = get_map(s)

  # Interpolate the map function to create a FEFunction
  Xₕ = interpolate_everywhere(map, U)

  return Xₕ
end

"""
  get_quasi_static_residual(s::MooringSegment, Xₕ::CellField, g::Real=9.81)

  Get the quasi-static residual for a segment. This function computes the residual based on the material properties,
  measures, and the reference configuration of the segment. It assumes that the segment is in a quasi-static state
  under the influence of gravity.
  Input:
  - `s::MooringSegment`: The segment for which the residual is computed.
  - `Xₕ::CellField`: The reference configuration of the segment.
  - `g::Real`: The gravitational acceleration (default is 9.81 m/s²).
  Output:
  - `res::Function`: The residual function that takes two arguments: the trial function `u` and the test function `v`.
"""
function get_quasi_static_residual(s::MooringSegment, Xₕ::CellField, g::Real=9.81)
  
  # Get the material properties
  material = get_material(s)

  # Get Measures
  dΩ, dΓ1, dΓ2 = get_measures(s)

  # TDC quantities
  J = TDC.J(Xₕ)
  Jabs = norm∘(J)
  Q = TDC.Q∘J

  # Define the stress functions
  S(u) = M.S(material, Xₕ, u)
  FΓ(u) = TDC.FΓ(u, Xₕ)
  K(u) = M.K∘(FΓ(u), S(u))

  # Define gravity force
  Fᵨ = VectorValue(0.0, -get_density(s) * get_area(s) * g)

  # Define the residual function
  res(u, v) = ∫(((∇(v)' ⋅ Q') ⊙ K(u)) * Jabs)dΩ -
              ∫((v ⋅ Fᵨ) * Jabs)dΩ

  return res
end

end 