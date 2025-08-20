module MooringTopology
using LinearAlgebra
using Gridap.TensorValues

export TopoPoint, TopoLine, MooringTopologyData, build_adjacency, assign_coords

"""
TopoPoint Struct

Topological point with coordinates. It includes the fields:
- id: Unique identifier for the point
- coords: Coordinates (x, y) or (x, y, z) of the point
- mesh_size: Size of the mesh element associated with the point
"""
struct TopoPoint
    id::Int
    tag::String
    coords::Vector{Float64}  # coordinates (x, y, z)  
    mesh_size::Float64
end

"""
TopoPoint(id::Int, coords::Vector{Float64}, mesh_size::Float64=1.0)

Constructor without tag and with default mesh size. The tag is defined as
"Point_[id]".
"""
function TopoPoint(id::Int, coords::Vector{Float64},mesh_size::Float64=1.0)
  tag = "Point_$id"
  return TopoPoint(id, tag, coords, mesh_size)
end

"""
get_id((p::TopoPoint)

Get the unique identifier of a topological point.
"""
get_id(p::TopoPoint) = p.id

"""
get_coords(p::TopoPoint)

Get the coordinates of a topological point.
"""
get_coords(p::TopoPoint) = p.coords

"""
TopoSegment Struct

Topological segment with start/end point ids and length. It includes the fields:
- id: Unique identifier for the segment
- start: ID of the start point
- stop: ID of the stop point  
- length: Length of the segment
"""
struct TopoSegment
    id::Int
    tag::String
    start::Int
    stop::Int
    length::Float64
end

"""
TopoSegment(id::Int, start::Int, stop::Int, length::Float64)

Constructor without tag. Default tag name is "Segment_[id]".
"""
function TopoSegment(id::Int, start::Int, stop::Int, length::Float64)
    tag = "Segment_$id"
    return TopoSegment(id, tag, start, stop, length)
end

"""
"""
get_id(s::TopoSegment) = s.id

"""
"""
get_start_point(s::TopoSegment) = s.start

"""
"""
get_stop_point(s::TopoSegment) = s.stop

"""
"""
get_length(s::TopoSegment) = s.length


"""
TopoloyData struct

This struct is used to store the topology data, including:
- points: Vector of topological points
- segments: Vector of topological segments
"""
struct MooringTopologyData
    points::Vector{TopoPoint}
    segments::Vector{TopoSegment}
end

"""
build_adjacency(topo::MooringTopologyData)

This function builds a graph representation of the topology, where each point
stores the IDs and lengths of all directly connected points. 

## Example
```julia
adj = build_adjacency(topo)
# Access neighbors of point with ID 1:
neighbors = adj[1]  # e.g. [(2, 10.0), (3, 5.5)]
```
"""
function build_adjacency(topo::MooringTopologyData)
    adj = Dict{Int, Vector{Tuple{Int, Float64}}}()
    for p in topo.points
        adj[p.id] = Vector{Tuple{Int, Float64}}()
    end
    for l in topo.segments
        push!(adj[l.start], (l.stop, l.length))
        push!(adj[l.stop], (l.start, l.length))
    end
    return adj
end

"""
assign_coords(topo::MooringTopologyData; anchor=1)

This function assigns 1D coordinates along the topological graph given the lengths of the segments.
It uses a depth-first search (DFS) approach to traverse the topology and assign coordinates.
"""
function assign_coords(topo::MooringTopologyData; anchor=1)
    adj = build_adjacency(topo)
    coords = Dict{Int, Float64}()
    visited = Set{Int}()
    
    coords[anchor] = 0.0
    stack = [anchor]
    
    while !isempty(stack)
        current = pop!(stack)
        push!(visited, current)
        
        for (nbr, length) in adj[current]
            if !haskey(coords, nbr)
                coords[nbr] = coords[current] + length
                push!(stack, nbr)
            end
        end
    end
    return coords
end

"""
get_physical_map(seg::TopoSegment, data::MooringTopologyData)

This function makes the physical map between two points of a segment.
Given a coordinate along the segment `r`, it returns the coordinate in the physical space.
"""
function get_physical_map(seg::TopoSegment, data::MooringTopologyData)
  p1 = data.points[get_start_point(seg)]
  p2 = data.points[get_stop_point(seg)]
  x_p1 = get_coords(p1)
  x_p2 = get_coords(p2)
  return function(r::VectorValue{1,Float64})
      t = r / norm(x_p2 .- x_p1)         # normalize parameter to [0,1]
      return VectorValue(x_p1 .+ t[1] .* (x_p2 .- x_p1)) # linear interpolation
  end
end

end # module
