module MooringLines
import Mooring.ParameterHandlers as PH
import Mooring.MooringSegments as Seg
import Mooring.MooringDiscreteModel as DM
import Mooring.PointMotion as PM
import Mooring.MooringPoints as Pts
import Mooring.Materials as Mat
using Gridap.TensorValues

export setup_lines

"""
MooringLine struct

This struct is used to define a mooring line in the mooring system. 
A mooring lines is defined by a set of segments, each segment is defined as [`MooringSegment`](@ref) types.
It includes the following fields:
- `segments::Vector{MooringSegment}`: Vector of segments that make up the mooring line
"""
struct MooringLine
  segments::Vector{Seg.MooringSegment}
end

"""
  get_segments(line::MooringLine)
  Get the segments of a mooring line.
"""
get_segments(line::MooringLine) = line.segments


function setup_lines(ph::PH.ParameterHandler)
  lines = MooringLine[]
  # for (line_id, line_params) in ph.lines
  #   segments = MooringSegment[]
  #   for seg_id in line_params.segments
  #     seg_params = ph.segments[seg_id]
  #     start_point = ph.points[seg_params.start_point]
  #     stop_point = ph.points[seg_params.stop_point]
  #     material = ph.materials[seg_params.material_tag]
  #     drag = haskey(ph.drags, seg_params.drag_tag) ? ph.drags[seg_params.drag_tag] : nothing
  #     seabed = haskey(ph.seabeds, seg_params.seabed_tag) ? ph.seabeds[seg_params.seabed_tag] : nothing
  #     segment = MooringSegment(seg_params.id, start_point, stop_point, material, seg_params.density,
  #                              seg_params.area, seg_params.length, drag, seabed)
  #     push!(segments, segment)
  #   end
  #   mooring_line = MooringLine(segments)
  #   push!(lines, mooring_line)
  # end
  # return lines

  # Extract line parameters
  ph_lines = ph.lines

  # Loop over lines
  for line in values(ph_lines)

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

  end
  return lines
end

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

end