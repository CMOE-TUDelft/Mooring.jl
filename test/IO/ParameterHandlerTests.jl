import Mooring.ParameterHandlers as PH

@testset "Default constructors" begin
  p = PH.PointParameters(1; coords=[1.0, 2.0])
  @test p.id == 1
  @test p.tag == "Point_1"
  @test p.coords == [1.0, 2.0]

  s = PH.SegmentParameters(1; start_point=1, stop_point=2)
  @test s.tag == "Segment_1"
  @test s.start_point == 1
  @test s.stop_point == 2

  l = PH.LineParameters(1; points=[1,2], segments=[1])
  @test l.tag == "Line_1"
  @test l.points == [1,2]

  d = PH.DragParameters(tag="drag1", dragType="Custom", nd=0.1, AStr=0.01)
  @test d.tag == "drag1"
  @test d.dragType == "Custom"
  @test d.nd == 0.1
  @test d.AStr == 0.01

  w = PH.WaveParameters(tag="sea_state", Hs=3.5, Tp=12.0)
  @test w.Hs == 3.5
  @test w.Tp == 12.0

  m = PH.MaterialParameters(tag="steel", type="LinearElastic", E=2e11, μ=0.3)
  @test m.tag == "steel"
  @test m.E == 2e11

  mo = PH.MotionParameters(tag="motion1", type="CustomMotion", f="sin(t)")
  @test mo.tag == "motion1"
  @test mo.f == "sin(t)"

  sb = PH.SeaBedParameters(tag="soft", kn=1000.0)
  @test sb.kn == 1000.0
end

@testset "ParameterHandler basics" begin
  ph = PH.ParameterHandler()
  
  # add some entries
  ph.points[1] = PH.PointParameters(1, coords=[0.0,0.0])
  ph.seabeds["soft"] = PH.SeaBedParameters(tag="soft", kn=1234.0)

  @test haskey(ph.points, 1)
  @test ph.points[1].tag == "Point_1"
  @test haskey(ph.seabeds, "soft")
  @test ph.seabeds["soft"].kn == 1234.0
end

@testset "Round-trip YAML (file)" begin
  ph = PH.ParameterHandler()
  ph.points[1] = PH.PointParameters(1, coords=[1.0,2.0])
  ph.seabeds["default_seabed"] = PH.SeaBedParameters()

  tmpfile = tempname() * ".yaml"
  PH.save_to_yaml(ph, tmpfile)
  ph2 = PH.load_from_yaml(tmpfile)

  @test ph2.points[1].coords == [1.0,2.0]
  @test haskey(ph2.seabeds, "default_seabed")
end

@testset "Round-trip JSON (file)" begin
  ph = PH.ParameterHandler()
  ph.materials["steel"] = PH.MaterialParameters(tag="steel", type="LinearElastic", E=2e11)

  tmpfile = tempname() * ".json"
  PH.save_to_json(ph, tmpfile)
  ph2 = PH.load_from_json(tmpfile)

  @test haskey(ph2.materials, "steel")
  @test ph2.materials["steel"].E == 2e11
end

@testset "Dict conversion" begin
  ph = PH.ParameterHandler()
  ph.points[1] = PH.PointParameters(1, coords=[3.0,4.0])
  dict = PH.ParameterHandler._handler_to_dict(ph)
  ph2 = PH.ParameterHandler._dict_to_handler(dict)
  @test ph2.points[1].coords == [3.0,4.0]
end

@testset "YAML in-memory round-trip" begin
  ph = PH.ParameterHandler()
  ph.waves["sea_state"] = PH.WaveParameters(tag="sea_state", Hs=5.0, Tp=15.0)

  dict = PH.ParameterHandler._handler_to_dict(ph)
  yaml_str = sprint(io -> YAML.write(io, dict))
  dict2 = YAML.load(yaml_str)
  ph2 = PH.ParameterHandler._dict_to_handler(dict2)

  @test haskey(ph2.waves, "sea_state")
  @test ph2.waves["sea_state"].Tp == 15.0
end

@testset "JSON in-memory round-trip" begin
  ph = PH.ParameterHandler()
  ph.drags["drag1"] = PH.DragParameters(tag="drag1", dragType="Custom", nd=0.05, AStr=0.002)

  dict = PH.ParameterHandler._handler_to_dict(ph)
  json_str = JSON3.write(dict; indent=2)
  dict2 = JSON3.read(json_str)
  ph2 = PH.ParameterHandler._dict_to_handler(dict2)

  @test haskey(ph2.drags, "drag1")
  @test ph2.drags["drag1"].nd == 0.05
end

@testset "Realistic example with nested YAML/JSON round-trip" begin
    ph = PH.ParameterHandler()

    # Define dependencies first
    ph.materials["steel"] = PH.MaterialParameters(tag="steel", type="LinearElastic", E=2e11, μ=0.3)
    ph.drags["drag1"] = PH.DragParameters(tag="drag1", dragType="Custom", nd=0.1, AStr=0.01, Cd_n=1.2)
    ph.seabeds["sand"] = PH.SeaBedParameters(tag="sand", kn=2000.0, linear_damping_factor=0.1)

    # Define points and segment linking those
    ph.points[1] = PH.PointParameters(1, coords=[0.0, 0.0])
    ph.points[2] = PH.PointParameters(2, coords=[10.0, 0.0])

    ph.segments[1] = PH.SegmentParameters(
        1,
        start_point=1,
        stop_point=2,
        length=10.0,
        density=7850.0,
        area=0.05,
        material_tag="steel",
        drag_tag="drag1",
        seabed_tag="sand"
    )

    ph.lines[1] = PH.LineParameters(1; points=[1,2], segments=[1])

    # --- YAML round-trip ---
    dict = PH.ParameterHandler._handler_to_dict(ph)
    yaml_str = sprint(io -> YAML.write(io, dict))
    dict2 = YAML.load(yaml_str)
    ph2 = PH.ParameterHandler._dict_to_handler(dict2)

    @test haskey(ph2.segments, 1)
    seg = ph2.segments[1]
    @test seg.material_tag == "steel"
    @test seg.drag_tag == "drag1"
    @test seg.seabed_tag == "sand"
    @test haskey(ph2.materials, "steel")
    @test haskey(ph2.drags, "drag1")
    @test haskey(ph2.seabeds, "sand")
    @test ph2.lines[1].segments == [1]

    # --- JSON round-trip ---
    json_str = JSON3.write(dict; indent=2)
    dict3 = JSON3.read(json_str)
    ph3 = PH.ParameterHandler._dict_to_handler(dict3)

    @test haskey(ph3.segments, 1)
    seg3 = ph3.segments[1]
    @test seg3.material_tag == "steel"
    @test ph3.materials["steel"].E == 2e11
    @test ph3.drags["drag1"].Cd_n == 1.2
    @test ph3.seabeds["sand"].kn == 2000.0
end