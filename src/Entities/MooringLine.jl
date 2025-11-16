module MooringLines
import Mooring.ParameterHandlers as PH
import Mooring.MooringSegments as Seg
import Mooring.MooringDiscreteModel as DM
import Mooring.PointMotion as PM
import Mooring.MooringPoints as Pts
import Mooring.Materials as Mat
using Gridap.TensorValues
using Gridap.FESpaces
using Gridap.MultiField
using Gridap.ODEs
using Gridap.Algebra
using Gridap.ReferenceFEs: num_point_dims
using Gridap.Geometry: get_node_coordinates
using Roots: find_zero

export setup_lines

"""
MooringLine struct

This struct is used to define a mooring line in the mooring system. 
A mooring lines is defined by a set of segments, each segment is defined as [`MooringSegment`](@ref) types.
It includes the following fields:
- `segments::Vector{MooringSegment}`: Vector of segments that make up the mooring line
"""
struct MooringLine
  segments::Dict{Int, Seg.MooringSegment}
end

"""
  get_segments(line::MooringLine)
  Get the segments of a mooring line.
"""
get_segments(line::MooringLine) = line.segments

"""
get_physical_map(seg::PH.SegmentParameters, ph::PH.ParameterHandler, start_point::Pts.MooringPoint, stop_point::Pts.MooringPoint)

This function makes the physical map between two points of a segment.
Given a coordinate along the segment `r`, it returns the coordinate in the physical space.
"""
function get_physical_map(seg::PH.SegmentParameters, ph::PH.ParameterHandler, start_point::Pts.MooringPoint, stop_point::Pts.MooringPoint)

  # Get physical coordinates of start and stop points
  p1_id = seg.start_point
  p2_id = seg.stop_point
  x_phys_p1 = ph.points[p1_id].coords
  x_phys_p2 = ph.points[p2_id].coords

  # Reference coordinates of start and stop points
  p1_ref_node = start_point.btrian.glue.face_to_bgface[1]
  p2_ref_node = stop_point.btrian.glue.face_to_bgface[1]
  x_ref_p1 = get_node_coordinates(start_point.btrian)[p1_ref_node][1]
  x_ref_p2 = get_node_coordinates(stop_point.btrian)[p2_ref_node][1]
  @assert norm(x_ref_p1-x_ref_p2) ≈ seg.length "Reference length between points $p1_id and $p2_id does not match segment length $(seg.length)"

  return function(r::VectorValue{1,Float64})
      s = (r[1] - x_ref_p1) / (x_ref_p2 - x_ref_p1)        # parametric coordinate along segment
      x_phys = x_phys_p1 .+ s .* (x_phys_p2 .- x_phys_p1) # linear interpolation
      return VectorValue(x_phys) 
  end
end

"""
get_physical_linear_map(seg::PH.SegmentParameters, ph::PH.ParameterHandler, start_point::Pts.MooringPoint, stop_point::Pts.MooringPoint)

This function creates a linear physical map between two points of a segment.
Given a reference coordinate along the segment `r` (from the 1D mesh), it returns the 
corresponding coordinate in 3D physical space using linear interpolation.

The mapping is constructed such that:
- The reference coordinates `x_ref_p1` and `x_ref_p2` from the 1D mesh are mapped to
- The physical coordinates `x_phys_p1` and `x_phys_p2` in 3D space
- The total reference length between points must match the segment's unstretched length
- The physical distance between endpoints can differ from the unstretched length (pre-tension/slack)

Arguments:
- `seg::PH.SegmentParameters`: Segment parameters including start/stop point IDs and length
- `ph::PH.ParameterHandler`: Parameter handler containing point coordinates
- `start_point::Pts.MooringPoint`: Start point with reference triangulation
- `stop_point::Pts.MooringPoint`: Stop point with reference triangulation

Returns:
- A function that maps `r::VectorValue{1,Float64}` (reference) to `VectorValue{D,Float64}` (physical)
"""
function get_physical_linear_map(seg::PH.SegmentParameters, ph::PH.ParameterHandler, start_point::Pts.MooringPoint, stop_point::Pts.MooringPoint)

  # Get physical coordinates of start and stop points
  p1_id = seg.start_point
  p2_id = seg.stop_point
  x_phys_p1 = ph.points[p1_id].coords
  x_phys_p2 = ph.points[p2_id].coords

  # Reference coordinates of start and stop points
  p1_ref_node = start_point.btrian.glue.face_to_bgface[1]
  p2_ref_node = stop_point.btrian.glue.face_to_bgface[1]
  x_ref_p1 = get_node_coordinates(start_point.btrian)[p1_ref_node][1]
  x_ref_p2 = get_node_coordinates(stop_point.btrian)[p2_ref_node][1]
  @assert norm(x_ref_p1-x_ref_p2) ≈ seg.length "Reference length between points $p1_id and $p2_id does not match segment length $(seg.length)"

  return function(r::VectorValue{1,Float64})
      s = (r[1] - x_ref_p1) / (x_ref_p2 - x_ref_p1)        # parametric coordinate s ∈ [0,1] along segment
      x_phys = x_phys_p1 .+ s .* (x_phys_p2 .- x_phys_p1) # linear interpolation in physical space
      return VectorValue(x_phys) 
  end
end

"""
get_physical_quadratic_map(seg::PH.SegmentParameters, ph::PH.ParameterHandler, start_point::Pts.MooringPoint, stop_point::Pts.MooringPoint)

This function creates a quadratic physical map between two points of a segment using a catenary-like curve.
Given a reference coordinate along the segment `r` (from the 1D mesh), it returns the 
corresponding coordinate in physical space.

The mapping is constructed by:
1. Computing the vertical plane containing the two endpoints
2. Solving for a parabolic catenary approximation in this 2D plane with arc length = seg.length
3. Rotating the curve back to 3D coordinates (or using directly in 2D)

For 2D problems, the vertical plane is simply the x-z plane.
For 3D problems, we find the vertical plane containing both points and compute the catenary there.

The parabola is defined as z(ξ) = -4*h*ξ*(1-ξ) where ξ ∈ [0,1] is the normalized horizontal coordinate
and h is the sag depth, computed such that the arc length equals seg.length.

Arguments:
- `seg::PH.SegmentParameters`: Segment parameters including start/stop point IDs and length
- `ph::PH.ParameterHandler`: Parameter handler containing point coordinates
- `start_point::Pts.MooringPoint`: Start point with reference triangulation
- `stop_point::Pts.MooringPoint`: Stop point with reference triangulation

Returns:
- A function that maps `r::VectorValue{1,Float64}` (reference) to `VectorValue{D,Float64}` (physical)
"""
function get_physical_quadratic_map(seg::PH.SegmentParameters, ph::PH.ParameterHandler, start_point::Pts.MooringPoint, stop_point::Pts.MooringPoint)

  # Get physical coordinates of start and stop points
  p1_id = seg.start_point
  p2_id = seg.stop_point
  x_phys_p1 = ph.points[p1_id].coords
  x_phys_p2 = ph.points[p2_id].coords

  # Reference coordinates of start and stop points
  p1_ref_node = start_point.btrian.glue.face_to_bgface[1]
  p2_ref_node = stop_point.btrian.glue.face_to_bgface[1]
  x_ref_p1 = get_node_coordinates(start_point.btrian)[p1_ref_node][1]
  x_ref_p2 = get_node_coordinates(stop_point.btrian)[p2_ref_node][1]
  @assert norm(x_ref_p1-x_ref_p2) ≈ seg.length "Reference length between points $p1_id and $p2_id does not match segment length $(seg.length)"

  # Determine dimensionality
  dims = length(x_phys_p1)
  
  # Target length is the segment length
  target_length = seg.length
  
  # Setup coordinate transformation to vertical plane
  if dims == 2
    # 2D case: already in x-z plane
    # Horizontal distance in plane
    horiz_dist_plane = abs(x_phys_p2[1] - x_phys_p1[1])
    vert_dist_plane = x_phys_p2[2] - x_phys_p1[2]  # Can be positive or negative
    
    # Transform function: 2D plane coords (ξ, z) -> physical coords (x, y)
    function plane_to_physical_2d(xi::Real, z_offset::Real)
      x = x_phys_p1[1] + xi * (x_phys_p2[1] - x_phys_p1[1])
      z = x_phys_p1[2] + xi * vert_dist_plane + z_offset
      return [x, z]
    end
    
    transform_to_physical = plane_to_physical_2d
    
  elseif dims == 3
    # 3D case: find vertical plane containing both points
    # The plane is defined by:
    # - Horizontal direction: from p1 to p2 projected onto x-y plane
    # - Vertical direction: z-axis
    
    # Horizontal vector in x-y plane
    horiz_vec_xy = [x_phys_p2[1] - x_phys_p1[1], x_phys_p2[2] - x_phys_p1[2], 0.0]
    horiz_dist_plane = norm(horiz_vec_xy)
    
    # Vertical distance
    vert_dist_plane = x_phys_p2[3] - x_phys_p1[3]
    
    # Unit vector in horizontal direction
    if horiz_dist_plane > 1e-10
      horiz_unit = horiz_vec_xy ./ horiz_dist_plane
    else
      # Points are vertically aligned, choose arbitrary horizontal direction
      horiz_unit = [1.0, 0.0, 0.0]
    end
    
    # Transform function: 2D plane coords (ξ, z_offset) -> 3D physical coords
    function plane_to_physical_3d(xi::Real, z_offset::Real)
      # Position along horizontal direction in the plane
      horiz_pos = xi * horiz_dist_plane
      # Vertical position (includes linear drop plus sag offset)
      vert_pos = x_phys_p1[3] + xi * vert_dist_plane + z_offset
      # 3D position
      x = x_phys_p1[1] + horiz_pos * horiz_unit[1]
      y = x_phys_p1[2] + horiz_pos * horiz_unit[2]
      z = vert_pos
      return [x, y, z]
    end
    
    transform_to_physical = plane_to_physical_3d
    
  else
    error("Unsupported dimensionality: $dims. Only 2D and 3D are supported.")
  end
  
  # Distance in the vertical plane
  plane_distance = sqrt(horiz_dist_plane^2 + vert_dist_plane^2)
  
  # If segment length equals plane distance, no sag needed
  if abs(target_length - plane_distance) < 1e-10
    sag_depth = 0.0
  else
    # Function to compute arc length of parabola z = -4*h*ξ*(1-ξ) for ξ ∈ [0,1]
    # in the vertical plane, combined with linear vertical drop
    function arc_length_in_plane(h::Real)
      # Use numerical integration (Simpson's rule)
      n = 1000
      dxi = 1.0 / n
      length_sum = 0.0
      
      for i in 0:n-1
        xi0 = i * dxi
        xi1 = (i + 1) * dxi
        xim = (xi0 + xi1) / 2.0
        
        # Position in plane at each point
        # Horizontal: ξ * horiz_dist_plane
        # Vertical: ξ * vert_dist_plane - 4*h*ξ*(1-ξ)
        
        # Derivatives
        # dhoriz/dξ = horiz_dist_plane (constant)
        dhoriz_dxi = horiz_dist_plane
        # dvert/dξ = vert_dist_plane - 4*h*(1-2*ξ)
        dvert_dxi_0 = vert_dist_plane - 4.0 * h * (1.0 - 2.0 * xi0)
        dvert_dxi_1 = vert_dist_plane - 4.0 * h * (1.0 - 2.0 * xi1)
        dvert_dxi_m = vert_dist_plane - 4.0 * h * (1.0 - 2.0 * xim)
        
        # Arc length element in plane: √((dhoriz/dξ)² + (dvert/dξ)²)
        dl_0 = sqrt(dhoriz_dxi^2 + dvert_dxi_0^2)
        dl_1 = sqrt(dhoriz_dxi^2 + dvert_dxi_1^2)
        dl_m = sqrt(dhoriz_dxi^2 + dvert_dxi_m^2)
        
        # Simpson's rule
        length_sum += (dl_0 + 4.0 * dl_m + dl_1) * dxi / 6.0
      end
      
      return length_sum
    end
    
    # Residual function: we want arc_length_in_plane(h) = target_length
    residual(h::Real) = arc_length_in_plane(h) - target_length
    
    # Initial guess for sag depth using catenary approximation
    excess_length = target_length - plane_distance
    h_initial = excess_length > 0 ? excess_length^2 / (8 * max(horiz_dist_plane, 1e-10)) : 0.0
    
    # Solve for sag depth using find_zero
    try
      sag_depth = find_zero(residual, h_initial)
    catch
      # If find_zero fails, use bisection with a reasonable bracket
      h_min = 0.0
      h_max = max(excess_length, horiz_dist_plane)
      try
        sag_depth = find_zero(residual, (h_min, h_max))
      catch
        # If still fails, no sag (use linear)
        @warn "Could not solve for sag depth, using linear interpolation"
        sag_depth = 0.0
      end
    end
  end

  return function(r::VectorValue{1,Float64})
      # Parametric coordinate along segment
      s = (r[1] - x_ref_p1) / (x_ref_p2 - x_ref_p1)  # s ∈ [0,1]
      
      # Compute position in vertical plane
      # Parabolic sag: z_offset = -4*h*s*(1-s)
      z_offset = -4.0 * sag_depth * s * (1.0 - s)
      
      # Transform from plane coordinates to physical coordinates
      x_result = transform_to_physical(s, z_offset)
      
      return VectorValue(x_result...) 
  end
end

"""
  setup_lines(ph::PH.ParameterHandler)

This function sets up the mooring lines based on the provided parameter handler.
It creates the discrete model, mooring points, and mooring segments for each line defined in the parameter handler.
It returns a dictionary of mooring lines, where the keys are the line IDs and the values are the corresponding `MooringLine` objects.
"""
function setup_lines(ph::PH.ParameterHandler)

  # Dictionary to store mooring lines
  lines = Dict{Int, MooringLine}()

  # Loop over lines
  for (line_id,line) in ph.lines

    # Create discrete model
    model = DM.generate_discrete_model(line, ph)

    # Create MooringPoints
    points = Dict{Int, Pts.MooringPoint}()
    for p_id in line.points
      point_params = ph.points[p_id]
      motion = PM.MotionType(p_id, ph)
      point = Pts.MooringPoint(model, point_params.tag, motion)
      points[p_id] = point
    end

    # Create MooringSegments
    segments = Dict{Int, Seg.MooringSegment}()
    for s_id in line.segments
      seg_params = ph.segments[s_id]
      start_point = points[seg_params.start_point]
      stop_point = points[seg_params.stop_point]
      map = get_physical_quadratic_map(seg_params, ph, start_point, stop_point)
      mat_params = ph.materials[seg_params.material_tag]
      material = Mat.Material(mat_params)
      seabed = ph.seabeds[seg_params.seabed_tag]
      segment = Seg.MooringSegment(model, 
                                   seg_params.tag,
                                   points[seg_params.start_point],
                                   points[seg_params.stop_point],
                                   map, 
                                   material, 
                                   seg_params.density,
                                   seg_params.area,
                                   seabed)
      segments[s_id] = segment
    end

    # Create MooringLine
    mooring_line = MooringLine(segments)

    # Store line
    lines[line_id] = mooring_line
  end

  return lines
end

"""
  get_Transient_FE_spaces(line::MooringLine)

This function retrieves the transient finite element spaces for all segments in a mooring line.
It returns a tuple containing:
- `X`: A `TransientMultiFieldFESpace` that combines the transient trial finite element spaces of all segments.
- `Y`: A `MultiFieldFESpace` that combines the test finite element spaces of all segments.
"""
function get_transient_FE_spaces(line::MooringLine)
  test_spaces = SingleFieldFESpace[]
  trial_spaces = TransientTrialFESpace[]
  for (s_id, segment) in line.segments
    dim = Seg.get_physical_dim(segment)
    U, V = Seg.get_transient_FESpaces(segment,dim=dim)
    push!(test_spaces, V)
    push!(trial_spaces, U)
  end
  X = TransientMultiFieldFESpace(trial_spaces)
  Y = MultiFieldFESpace(test_spaces)
  return X, Y
end

"""
  get_reference_configuration(line::MooringLine, X::MultiFieldFESpace)

This function computes the reference configuration for a mooring line given a multi-field finite element space `X`.
It returns a finite element function `Xₕ` that represents the reference configuration of the mooring line.
The reference configuration is obtained by interpolating the physical maps of each segment in the mooring line
over the provided finite element space."""
function get_reference_configuration(line::MooringLine, X::MultiFieldFESpace)
  maps = Function[]
  for (s_id, segment) in line.segments
    map = Seg.get_map(segment)
    push!(maps, map)
  end
  Xₕ = interpolate_everywhere(maps, X)
  return Xₕ
end

"""
  get_quasi_static_residual(line::MooringLine, Xₕ::MultiFieldFEFunction; g::Float64=9.81)

This function computes the quasi-static residual for a mooring line given a multi-field finite element function `Xₕ`.
It returns a function `res` that represents the quasi-static residual of the mooring line.
The residual is computed by summing the quasi-static residuals of each segment in the mooring line.
An optional gravitational acceleration parameter `g` can be provided (default is 9.81 m/s²).

Check also [`MooringSegment.get_quasi_static_residual`](@ref) for details on the definition of the segment residual.
"""
function get_quasi_static_residual(line::MooringLine, Xₕ::MultiFieldFEFunction, g::Float64=9.81)
  res_terms = Function[]
  for (s_it,(s_id, segment)) in enumerate(line.segments)
    res = Seg.get_quasi_static_residual(segment, Xₕ[s_it], g)
    push!(res_terms, res)
  end
  function res(x,y)
    sum(res_terms[i](x[i],y[i]) for i in eachindex(res_terms))
  end
  return res
end

"""
  solve_quasistatic(ph::PH.ParameterHandler)

This function solves the quasi-static problem for all mooring lines defined in the provided parameter handler `ph`.
It sets up the mooring lines, computes the transient finite element spaces, reference configurations, and
quasi-static residuals for each line, and then solves the resulting nonlinear system using a Newton solver.

It returns a vector of finite element functions `x`, where each function corresponds to the solution of a mooring line.

Check also [`get_transient_FE_spaces`](@ref), [`get_reference_configuration`](@ref), and [`get_quasi_static_residual`](@ref) for details on the individual steps.
"""
function solve_quasistatic(ph::PH.ParameterHandler)
  # Setup lines
  mlines = setup_lines(ph)

  u = FEFunction[]
  x_ref = FEFunction[]

  # Loop over lines
  for (line_id,line) in mlines

    @info("Solving line $line_id with $(length(line.segments)) segments")

    # Transient FE spaces
    X, Y = get_transient_FE_spaces(line)

    # Reference configuration
    Xₕ = get_reference_configuration(line, X(0.0))

    # Quasi-static residual
    # TODO: add g as input parameter (global constants)
    res = get_quasi_static_residual(line, Xₕ)
  
    # solve
    # TODO: add nls parameters as input parameters (solver parameters)
    op = FEOperator(res, X(0.0), Y)
    nls = NLSolver(BackslashSolver(), iterations=200, show_trace=true, ftol=1e-8,method=:newton)
    uₕ = solve(nls, op)
    push!(u, uₕ)
    push!(x_ref, Xₕ)

  end

  return u,x_ref
end

end # module