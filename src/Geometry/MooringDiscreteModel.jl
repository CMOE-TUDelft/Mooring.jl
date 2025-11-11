module MooringDiscreteModel

using GridapGmsh: gmsh, GmshDiscreteModel
import Mooring.ParameterHandlers as PH

"""
build_adjacency(line::PH.LineParameters, ph::PH.ParameterHandler)

This function builds a graph representation of the line topology, where each point
stores the IDs and lengths of all directly connected points. 

## Example
```julia
adj = build_adjacency(line, ph)
# Access neighbors of point with ID 1:
neighbors = adj[1]  # e.g. [(2, 10.0), (3, 5.5)]
```
"""
function build_adjacency(line::PH.LineParameters, ph::PH.ParameterHandler)
    adj = Dict{Int, Vector{Tuple{Int, Float64}}}()
    for p_id in line.points
        adj[p_id] = Vector{Tuple{Int, Float64}}()
    end
    for s_id in line.segments
        s = ph.segments[s_id]
        push!(adj[s.start_point], (s.stop_point, s.length))
        push!(adj[s.stop_point], (s.start_point, s.length))
    end
    return adj
end

"""
assign_coords(line::PH.LineParameters, ph::PH.ParameterHandler; anchor=1)

This function assigns 1D coordinates along the topological graph given the lengths of the segments.
It uses a depth-first search (DFS) approach to traverse the topology and assign coordinates.
"""
function assign_coords(line::PH.LineParameters, ph::PH.ParameterHandler; anchor=1)
    adj = build_adjacency(line, ph)
    coords = Dict{Int, Float64}()
    visited = Set{Int}()

    # Check anchor exists
    if !haskey(adj, anchor)
        println("Anchor point $anchor not found in line points. Assigning first point as anchor.")
        anchor = first(keys(adj))
    end

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
generate_mesh(line::PH.LineParameters,ph::PH.ParameterHandler)

Generate a 1D mesh of the topological graph defined using the ParameterHandler. 
The mesh size can be controlled by the `mesh_size` parameter in the PointParameters.
The function returns a `gmsh` model. 
"""
function generate_mesh(line::PH.LineParameters, ph::PH.ParameterHandler)
  gmsh.initialize()
  gmsh.model.add("mooring_model")

  # 1D topological graph coordinates
  coords_dict = assign_coords(line, ph)

  # Add points
  for p_id in line.points
    p = ph.points[p_id]
    coords = coords_dict[p_id]
    gmsh.model.geo.addPoint(coords, 0.0, 0.0, p.mesh_size, p_id)
    gmsh.model.addPhysicalGroup(0, [p_id], p_id)
    gmsh.model.setPhysicalName(0, p_id, p.tag)
  end

  # Add lines
  for s_id in line.segments
    s = ph.segments[s_id]
    gmsh.model.geo.addLine(s.start_point, s.stop_point, s_id)
    gmsh.model.addPhysicalGroup(1, [s_id], s_id)
    gmsh.model.setPhysicalName(1, s_id, s.tag)
  end

  gmsh.model.geo.synchronize()
  gmsh.model.mesh.generate(1)

  return gmsh
end

"""
generate_discrete_model(line::PH.LineParameters, ph::PH.ParameterHandler)

This function generates a discrete model from the given line parameters and parameter handler.
"""
function generate_discrete_model(line::PH.LineParameters, ph::PH.ParameterHandler)
  gmsh = generate_mesh(line, ph)
  model = GmshDiscreteModel(gmsh)
  gmsh.finalize()
  return model
end

end # module MooringDiscreteModel