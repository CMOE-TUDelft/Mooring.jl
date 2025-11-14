module ThreeLine
import Mooring.ParameterHandlers as PH
import Mooring.MooringLines as ML
using Gridap.Geometry
using Gridap.Visualization: writevtk
using LinearAlgebra

# System geometry parameters (following MoorPy example)
depth = 600.0                                   # water depth [m]
angles = [π/3, π, 5π/3]                         # line headings [rad] (60°, 180°, 300°)
rAnchor = 1600.0                                # anchor radius/spacing [m]
zFair = -21.0                                   # fairlead z elevation [m]
rFair = 20.0                                    # fairlead radius [m]
lineLength = 2100.0                             # line unstretched length [m]

# Create ParameterHandler
ph = PH.ParameterHandler()

# Define material properties (120mm chain equivalent)
ph.materials["chain1"] = PH.MaterialParameters(
  tag="chain1", 
  E=2.1e11,      # Steel Young's modulus [Pa]
  μ=8.1e10       # Steel shear modulus [Pa]
)

# Define motion types
ph.motions["fixed"] = PH.MotionParameters(tag="fixed",
  type="CustomMotion", 
  f="(t,x)->VectorValue(0.0, 0.0, 0.0)"  # Initially at equilibrium
)
ph.motions["floating_body"] = PH.MotionParameters(
  tag="floating_body", 
  type="CustomMotion", 
  f="(t,x)->VectorValue(0.0, 0.0, 0.0)"  # Initially at equilibrium
)

# Loop over lines
for (i,angle) in enumerate(angles)
  
  # Create anchor and fairlead points for each line
  ph.points[2*(i-1)+1] = PH.PointParameters(
    id=2*(i-1)+1, 
    coords=[rAnchor*cos(angle), rAnchor*sin(angle), -depth], 
    motion_tag="fixed", 
    mesh_size=50.0
  )
  ph.points[2i] = PH.PointParameters(
    id=2i, 
    coords=[rFair*cos(angle), rFair*sin(angle), zFair], 
    motion_tag="floating_body", 
    mesh_size=50.0
  )
  
  # Define segment for each line
  ph.segments[i] = PH.SegmentParameters(
    id=i, 
    start_point=2*(i-1)+1, 
    stop_point=2i, 
    material_tag="chain1", 
    area=π*(0.06)^2,          # Approximate area for 120mm chain [m²]
    density=7850.0,           # Steel density [kg/m³]
    length=lineLength, 
    drag_tag="no_drag", 
    seabed_tag="default_seabed"
  )
  
  # Define line (three un-connected lines)
  ph.lines[i] = PH.LineParameters(
    id=i, 
    points=[2*(i-1)+1, 2i], 
    segments=[i]
  )
  
end

# Solve quasi-static equilibrium
u₀,x_ref = ML.solve_quasistatic(ph)

for (iline, (u_line, x_ref_line)) in enumerate(zip(u₀, x_ref))
  writevtk(get_triangulation(u_line), 
           "three_line_deformed_line$(iline).vtu", 
           cellfields=["Displacement"=>u_line[1], "ReferenceConfiguration"=>x_ref_line[1]])
end


end