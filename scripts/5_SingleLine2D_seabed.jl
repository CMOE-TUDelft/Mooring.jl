module SingleLine2D
import Mooring.ParameterHandlers as PH
import Mooring.MooringLines as ML
using Gridap.Geometry
using Gridap.Visualization: writevtk

# System geometry parameters
depth = 0.5                   # water depth [m]
lineLength = 1.1540019835742  # line unstretched length [m]
anchor_x = 1.0                # anchor x coordinate [m]
anchor_y = 0.0                # anchor y coordinate [m]
fairlead_y = 0.5              # fairlead y elevation [m]
fairlead_x = 0.0              # fairlead x coordinate [m]

# Create ParameterHandler
ph = PH.ParameterHandler()

# Define material properties
ph.materials["line_material"] = PH.MaterialParameters(
  tag="line_material",
  E=1.0e4,        # Young's modulus [Pa]
  μ=5.0e3,         # Shear modulus [Pa]
)

# Define motion types
ph.motions["fixed"] = PH.MotionParameters(
  tag="fixed",
  type="CustomMotion",
  f="(t,x)->VectorValue(0.0, 0.0)"
)  # Fixed point

# Define points
ph.points[1] = PH.PointParameters(
  id=1,
  coords = [fairlead_x, fairlead_y],
  motion_tag="fixed",
  mesh_size=0.01
)
ph.points[2] = PH.PointParameters(
  id=2,
  coords = [anchor_x, anchor_y],
  motion_tag="fixed",
  mesh_size=0.01
)

# Define seabed parameters
area = 0.01
ρ = 2000.0/(9.81*area)
weight_per_length = ρ * 9.81 * area
ph.seabeds["default_seabed"] = PH.SeaBedParameters(
  tag="default_seabed",
  kn=10000.0,    # seabed stiffness [N/m]
  still_weight = weight_per_length
)

# Define segment
ph.segments[1] = PH.SegmentParameters(
  id=1,
  tag="line1",
  start_point=1,
  stop_point=2,
  material_tag="line_material",
  area=0.01,            # Cross-sectional area [m²]
  density=ρ,
  length=lineLength,
  seabed_tag="default_seabed"
)

# Define line
ph.lines[1] = PH.LineParameters(
  id=1,
  points=[1,2],
  segments=[1]
)

# Solve quasi-static problem
u, x = ML.solve_quasistatic(ph)

for (iline, (u_line, x_ref_line)) in enumerate(zip(u,x))
  writevtk(get_triangulation(u_line), 
           "2d_line_$(iline).vtu", 
           cellfields=["u"=>u_line[1], "x"=>x_ref_line[1]])
end

end # module SingleLine2D