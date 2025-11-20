module SingleLineMultiSegment2D
import Mooring.ParameterHandlers as PH
import Mooring.MooringLines as ML
using Gridap.Geometry
using Gridap.Visualization: writevtk

# System geometry parameters
depth = 0.5                   # water depth [m]
anchor_x = 1.0                # anchor x coordinate [m]
anchor_y = 0.0                # anchor y coordinate [m]
fairlead_y = 0.5              # fairlead y elevation [m]
fairlead_x = 0.0              # fairlead x coordinate [m]

# Segment lengths
chain1_length = 0.4           # Bottom chain segment [m]
polyester_length = 0.55        # Middle polyester segment [m]
chain2_length = 0.4           # Top chain segment [m]
total_length = chain1_length + polyester_length + chain2_length

# Create ParameterHandler
ph = PH.ParameterHandler()

# Define material properties
# Chain material (steel-like, stiffer)
ph.materials["chain"] = PH.MaterialParameters(
  tag="chain",
  E=2.1e11,       # Steel Young's modulus [Pa]
  μ=8.1e10        # Steel shear modulus [Pa]
)

# Polyester material (more flexible)
ph.materials["polyester"] = PH.MaterialParameters(
  tag="polyester",
  E=5.0e9,        # Polyester Young's modulus [Pa]
  μ=2.0e9         # Polyester shear modulus [Pa]
)

# Define motion types
ph.motions["fixed"] = PH.MotionParameters(
  tag="fixed",
  type="CustomMotion",
  f="(t,x)->VectorValue(0.0, 0.0)"
)  # Fixed point
ph.motions["free"] = PH.MotionParameters(
  tag="free",
  type=nothing
)  # Free point

# Define points (4 points: anchor, 2 intermediate, fairlead)
ph.points[1] = PH.PointParameters(
  id=1,
  coords = [fairlead_x, fairlead_y],
  motion_tag="fixed",
  mesh_size=0.1
)
ph.points[2] = PH.PointParameters(
  id=2,
  coords = [fairlead_x + 0.3, fairlead_y - 0.1],  # Intermediate point 1
  motion_tag="free",
  mesh_size=0.1
)
ph.points[3] = PH.PointParameters(
  id=3,
  coords = [fairlead_x + 0.7, anchor_y + 0.1],    # Intermediate point 2
  motion_tag="free",
  mesh_size=0.1
)
ph.points[4] = PH.PointParameters(
  id=4,
  coords = [anchor_x, anchor_y],
  motion_tag="fixed",
  mesh_size=0.1
)

# Define segments with different materials
# Segment 1: Top chain (from fairlead to first intermediate)
chain_area = 0.001            # Chain cross-sectional area [m²]
chain_density = 7850.0        # Steel density [kg/m³]
ph.segments[1] = PH.SegmentParameters(
  id=1,
  tag="chain_top",
  start_point=1,
  stop_point=2,
  material_tag="chain",
  area=chain_area,
  density=chain_density,
  length=chain2_length
)

# Segment 2: Middle polyester (between intermediate points)
polyester_area = 0.002        # Polyester cross-sectional area [m²]
polyester_density = 1380.0    # Polyester density [kg/m³]
ph.segments[2] = PH.SegmentParameters(
  id=2,
  tag="polyester",
  start_point=2,
  stop_point=3,
  material_tag="polyester",
  area=polyester_area,
  density=polyester_density,
  length=polyester_length
)

# Segment 3: Bottom chain (from second intermediate to anchor)
ph.segments[3] = PH.SegmentParameters(
  id=3,
  tag="chain_bottom",
  start_point=3,
  stop_point=4,
  material_tag="chain",
  area=chain_area,
  density=chain_density,
  length=chain1_length
)

# Define line with three segments
ph.lines[1] = PH.LineParameters(
  id=1,
  points=[1,2,3,4],
  segments=[1,2,3]
)

# Solve quasi-static problem
u, x = ML.solve_quasistatic(ph)

for (iline, (u_line, x_ref_line)) in enumerate(zip(u,x))
  for i_segment in 1:length(u_line)
    Ω = get_triangulation(u_line[i_segment])
    writevtk(Ω, "multisegment_2d_line_$(iline)_segment_$(i_segment).vtu", 
             cellfields=["u"=>u_line[i_segment], "x"=>x_ref_line[i_segment]])
  end
end

end # module SingleLineMultiSegment2D