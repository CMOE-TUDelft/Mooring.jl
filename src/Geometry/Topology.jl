module Topology

export TopoPoint, TopoLine, TopologyData, build_adjacency, assign_coords

"""
TopoPoint Struct

Topological point with coordinates. It includes the fields:
- id: Unique identifier for the point
- coords: Coordinates (x, y) or (x, y, z) of the point
"""
struct TopoPoint
    id::Int
    coords::Vector{Float64}  # coordinates (x, y, z)  
end

"""
TopoSegment Struct

Topological segment with start/end point ids, length, and physical mapping. It includes the fields:
- id: Unique identifier for the segment
- start: ID of the start point
- stop: ID of the stop point  
- length: Length of the segment
- physmap: Physical mapping function
"""
struct TopoSegment
    id::Int
    start::Int
    stop::Int
    length::Float64
    physmap::Function   # function r in [0,length] -> (x,y,z)
end

"""
TopoloyData struct

This struct is used to store the topology data, including:
- points: Vector of topological points
- segments: Vector of topological segments
"""
struct TopologyData
    points::Vector{TopoPoint}
    segments::Vector{TopoSegment}
end

"""
build_adjacency(topo::TopologyData)

This function builds a graph representation of the topology, where each point
stores the IDs and lengths of all directly connected points. 

## Example
```julia
adj = build_adjacency(topo)
# Access neighbors of point with ID 1:
neighbors = adj[1]  # e.g. [(2, 10.0), (3, 5.5)]
```
"""
function build_adjacency(topo::TopologyData)
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
assign_coords(topo::TopologyData; anchor=1)

This function assigns 1D coordinates along the topological graph given the lengths of the segments.
It uses a depth-first search (DFS) approach to traverse the topology and assign coordinates.
"""
function assign_coords(topo::TopologyData; anchor=1)
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

end # module
