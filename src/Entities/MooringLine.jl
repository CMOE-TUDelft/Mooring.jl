module MooringLines

"""
MooringLine struct

This struct is used to define a mooring line in the mooring system. 
A mooring lines is defined by a set of segments, each segment is defined as [`MooringSegment`](@ref) types.
It includes the following fields:
- `segments::Vector{MooringSegment}`: Vector of segments that make up the mooring line
"""
struct MooringLine
  segments::Vector{MooringSegment}
end

"""
  get_segments(line::MooringLine)
  Get the segments of a mooring line.
"""
get_segments(line::MooringLine) = line.segments


end