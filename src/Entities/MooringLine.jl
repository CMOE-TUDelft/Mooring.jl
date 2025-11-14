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
      s = (r[1] -x_ref_p1) / (x_ref_p2 - x_ref_p1)        # parametric coordinate along segment
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

This function creates a quadratic physical map between two points of a segment.
Given a reference coordinate along the segment `r` (from the 1D mesh), it returns the 
corresponding coordinate in 3D physical space using quadratic interpolation with a catenary-like shape.

The mapping is constructed such that:
- The reference coordinates `x_ref_p1` and `x_ref_p2` from the 1D mesh are mapped to
- The physical coordinates `x_phys_p1` and `x_phys_p2` in 3D space
- A parabolic sag is added in the vertical direction to approximate a catenary shape
- The sag depth is proportional to the horizontal span and segment length
- The total reference length between points must match the segment's unstretched length

The quadratic map provides a better initial guess for hanging cables compared to linear interpolation,
reducing the number of Newton iterations needed for convergence.

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

  # Calculate horizontal span and sag parameters
  dims = length(x_phys_p1)
  horizontal_span = norm(x_phys_p2[1:dims-1] .- x_phys_p1[1:dims-1])
  
  # Estimate sag depth based on cable length excess over horizontal span
  # For a catenary, the sag is approximately (L²-H²)/(8H) where L is cable length and H is horizontal span
  excess_length = seg.length - horizontal_span
  sag_depth = excess_length > 0 ? excess_length^2 / (8 * max(horizontal_span, 1.0)) : 0.0

  return function(r::VectorValue{1,Float64})
      s = (r[1] - x_ref_p1) / (x_ref_p2 - x_ref_p1)        # parametric coordinate s ∈ [0,1] along segment
      
      # Linear interpolation in horizontal directions
      x_phys = x_phys_p1 .+ s .* (x_phys_p2 .- x_phys_p1)
      
      # Add parabolic sag in vertical direction (last dimension)
      # Parabola: y = -4*h*s*(1-s) where h is the sag depth
      # Maximum sag occurs at s=0.5
      vertical_offset = -4.0 * sag_depth * s * (1.0 - s)
      
      # Create result vector with sag applied to vertical component
      x_result = collect(x_phys)
      x_result[dims] += vertical_offset
      
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
      segment = Seg.MooringSegment(model, 
                                   seg_params.tag,
                                   points[seg_params.start_point],
                                   points[seg_params.stop_point],
                                   map, 
                                   material, 
                                   seg_params.density,
                                   seg_params.area)
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