module MooringDiscreteModel

using GridapGmsh: gmsh, GmshDiscreteModel
import Mooring.ParameterHandlers as PH

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
  coords_dict = Topo.assign_coords(line, ph)

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