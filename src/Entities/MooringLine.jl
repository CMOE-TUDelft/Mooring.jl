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
get_physical_map(seg::PH.SegmentParameters, ph::PH.ParameterHandler)

This function makes the physical map between two points of a segment.
Given a coordinate along the segment `r`, it returns the coordinate in the physical space.
"""
function get_physical_map(seg::PH.SegmentParameters, ph::PH.ParameterHandler)
  p1 = seg.start_point
  p2 = seg.stop_point
  x_p1 = ph.points[p1].coords
  x_p2 = ph.points[p2].coords
  return function(r::VectorValue{1,Float64})
      t = r / norm(x_p2 .- x_p1)         # normalize parameter to [0,1]
      return VectorValue(x_p1 .+ t[1] .* (x_p2 .- x_p1)) # linear interpolation
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
      map = get_physical_map(seg_params, ph)
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
    U, V = Seg.get_transient_FESpaces(segment)
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

function solve_quasistatic(ph::PH.ParameterHandler)
  # Setup lines
  mlines = setup_lines(ph)

  x = FEFunction[]

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
    nls = NLSolver(BackslashSolver(), iterations=200)
    xₕ = solve(nls, op)
    push!(x, xₕ)

  end

  return x
end

end # module