module MooringDiscreteModel

using GridapGmsh: gmsh, GmshDiscreteModel
import Mooring.Topology as Topo

"""
generate_mesh(topo::Topo.TopologyData)

Generate a 1D mesh of the topological graph defined in `topo`. The mesh size
can be controlled by the `mesh_size` parameter in the `TopoPoint` struct.
The function returns a `gmsh` model. 
"""
function generate_mesh(topo::Topo.TopologyData)   
  gmsh.initialize()
  gmsh.model.add("mooring_model")

  # 1D topological graph coordinates
  coords_dict = Topo.assign_coords(topo)

  # Add points
  for p in topo.points
    coords = coords_dict[p.id]
    gmsh.model.geo.addPoint(coords, 0.0, 0.0, p.mesh_size, p.id)
    gmsh.model.addPhysicalGroup(0, [p.id], p.id)
    gmsh.model.setPhysicalName(0, p.id, p.tag)
  end

  # Add lines
  for l in topo.segments
      gmsh.model.geo.addLine(l.start, l.stop, l.id)
      gmsh.model.addPhysicalGroup(1, [l.id], l.id)
      gmsh.model.setPhysicalName(1, l.id, l.tag)
  end

  gmsh.model.geo.synchronize()
  gmsh.model.mesh.generate(1)

  return gmsh
end

"""
TO DO -> generate_discrete_model
"""

function generate_discrete_model(topo::Topo.TopologyData)
  gmsh = generate_mesh(topo)
  model = GmshDiscreteModel(gmsh)
  gmsh.finalize()
  return model
end

end # module MooringDiscreteModel