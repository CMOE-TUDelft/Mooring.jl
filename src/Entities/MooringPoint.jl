module MooringPoints
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.CellData: lazy_map, IntegrationMap, CellFieldAt
using Gridap.Geometry: get_faces, get_node_coordinates, CompositeTriangulation
using Gridap.TensorValues: VectorValue
using Mooring.PointMotion: MotionType, get_point_motion_function
import Mooring.TangentialDiffCalculus as TDC
import Mooring.Materials as M

MooringPointMotion = Union{MotionType, Nothing}

"""
MooringPoint Struct
This struct is used to define a point in the mooring system.
It includes the following fields:
- `tag::String`: Identifier for the point
- `btrian::Triangulation`: BoundaryTriangulation or InterfaceTriangulation with the position of the point in space
- `s_id_trians::Vector{Tuple{Int,Triangulation}}`: Vector of tuples containing segment ID and corresponding Triangulation connected to this point
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
  get_segment_ids(p::MooringPoint)
  Get the segment IDs connected to a point (ordered as in s_id_trians)
"""
get_segment_ids(p::MooringPoint) = [s_id_trian[1] for s_id_trian in p.s_id_trians]

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

is_free_point(p::MooringPoint) = p.motion === nothing

"""
get_measure(p::MooringPoint, degree::Int=1)
Get the measure for a mooring point. This function returns the integration measure for the point's triangulation.
In practice this is equivalent to evaluating the function at the point (0D entity).
"""
get_measure(p::MooringPoint, degree::Int=1) = Measure(get_triangulation(p), degree)

"""
  get_quasi_static_residual(p::MooringPoint, Xₕ::Vector{CellField}, g::Float64=9.81)
  Get the quasi-static residual contribution of a mooring point.
  This function computes the residual contribution of the point based on its motion type
  and the reference configuration provided by the vector of CellField `Xₕ`, which
  should correspond to the segments connected to the point.
  Input:
  - `p::MooringPoint`: The mooring point for which the residual is computed.
  - `materials::Vector{M.Material}`: Vector of Material corresponding to the segments connected to the point.
  - `Xₕ::Vector{CellField}`: Vector of CellField representing the reference configurations of the segments connected to the point.
  - `g::Float64`: Gravity acceleration (default is 9.81 m/s²).
Output:
  - `residual_function::Function`: A function that computes the residual contribution of the point.
"""
function get_quasi_static_residual(p::MooringPoint, materials::Vector{M.Material}, Xₕ::Vector{<:SingleFieldFEFunction}, g::Float64=9.81)
  
  # Get triangulation and measure (integration point)
  dΓ = get_measure(p)
  Γ = get_triangulation(p)
  points = get_cell_points(dΓ.quad)
  weights = dΓ.quad.cell_weight
  
  # TDC quantities at the point
  Xh1 = CellFieldAt{:plus}(Xₕ[1])
  Xh2 = CellFieldAt{:minus}(Xₕ[2])
  T1 = TDC.T∘(TDC.J(Xh1))
  T2 = TDC.T∘(TDC.J(Xh2))

  # Residual function
  # res = ∫([𝒘⋅ 𝐍𝜕Γ𝑿] ⋅ {𝐊})dΓ + ∫({𝒘}⋅ [𝐊 ⋅ 𝐍𝜕Γ𝑿])dΓ + ∫([u⋅ 𝐍𝜕Γ𝑿] ⋅ {𝐊})dΓ + ∫([𝒘] ⋅ [u])dΓ 
  # res = ∫([𝒘⋅ 𝐍𝜕Γ𝑿] ⋅ {𝐊})dΓ + ∫({𝒘}⋅     0      )dΓ + ∫([u⋅ 𝐍𝜕Γ𝑿] ⋅ {𝐊})dΓ + ∫([𝒘] ⋅ [u])dΓ 
  function res((u_l,u_r),(v_l,v_r)) 

    # Auxiliary quantities
    jump_u = u_l.⁺ - u_r.⁻
    jump_v = v_l.⁺ - v_r.⁻
    mean_v = (v_l.⁺ + v_r.⁻)/2

    # Force contributions from side 1
    FΓ1 = TDC.FΓ(u_l.⁺, Xh1)
    S1 = M.S(materials[1], Xh1, u_l.⁺)
    K1 = CellFieldAt{:plus}(M.K∘(FΓ1, S1))
    H1 = K1⋅T1

    # Force contributions from side 2
    FΓ2 = TDC.FΓ(u_r.⁻, Xh2)
    S2 = M.S(materials[2], Xh2, u_r.⁻)
    K2 = CellFieldAt{:minus}(M.K∘(FΓ2, S2))
    H2 = K2⋅T2
    
    # Domain contribution (TODO: hard-coded constant 10000 for penalty)
    c = DomainContribution()
    add_contribution!(c,Γ, lazy_map(IntegrationMap(),(10000*(jump_u⋅jump_v))(points), weights))
    add_contribution!(c,Γ, lazy_map(IntegrationMap(),((mean_v⋅(H1-H2))⋅VectorValue(1.0))(points), weights))
    return c

  end

  return res
end

end