# module ParameterHandler

# using YAML
# using JSON

# export ParameterHandler, from_defaults, from_yaml, from_json, from_dat_folder,
#        to_yaml, to_json, get_parameter, set_parameter

# # --- Central Container ---
# struct ParameterHandler
#     points::Dict{Int, Point}
#     segments::Dict{Int, Segment}
#     lines::Dict{Int, Line}
#     drag::Dict{String, DragProperties}
#     waves::WaveParameters
#     materials::Dict{String, Material}
#     motions::Dict{Int, Motion}
#     seabed::SeaBedParameters
# end

# # --- Constructors ---

# """
#     from_defaults() -> ParameterHandler

# Initialize a ParameterHandler with default values.
# """
# function from_defaults()
#     return ParameterHandler(
#         Dict{Int, Point}(),
#         Dict{Int, Segment}(),
#         Dict{Int, Line}(),
#         Dict{String, DragProperties}(),
#         WaveParameters(),  # default constructor
#         Dict{String, Material}(),
#         Dict{Int, Motion}(),
#         SeaBedParameters()
#     )
# end

# """
#     from_yaml(path::String) -> ParameterHandler

# Parse YAML input into a ParameterHandler.
# """
# function from_yaml(path::String)
#     data = YAML.load_file(path)
#     return parse_dict(data)
# end

# """
#     from_json(path::String) -> ParameterHandler
# """
# function from_json(path::String)
#     data = JSON.parsefile(path)
#     return parse_dict(data)
# end

# """
#     from_dat_folder(folder::String) -> ParameterHandler

# Aggregate multiple .dat files into a single ParameterHandler and 
# produce a YAML export.
# """
# function from_dat_folder(folder::String)
#     # (placeholder) implement file aggregation
#     data = Dict()  # aggregated dictionary from .dat files
#     handler = parse_dict(data)

#     # export aggregated version for reproducibility
#     to_yaml(handler, joinpath(folder, "aggregated.yaml"))
#     return handler
# end

# # --- Export ---
# function to_yaml(ph::ParameterHandler, path::String)
#     open(path, "w") do io
#         YAML.write(io, to_dict(ph))
#     end
# end

# function to_json(ph::ParameterHandler, path::String)
#     open(path, "w") do io
#         JSON.print(io, to_dict(ph))
#     end
# end

# # --- Internal Conversion Utilities ---
# function parse_dict(data::Dict)::ParameterHandler
#     # Map Dict to strongly typed structs
#     points = Dict(id => Point(id, d["tag"], d["coordinates"], d["motion"], d["mesh_size"]) for (id, d) in data["points"])
#     segments = Dict(id => Segment(id, d["tag"], d["start"], d["stop"], d["length"],
#                                   d["material"], d["density"], d["area"]) for (id, d) in data["segments"])
#     lines = Dict(id => Line(d["points"], d["segments"]) for (id, d) in data["lines"])
#     drag = Dict(name => DragProperties(; d... ) for (name, d) in data["drag"])
#     waves = WaveParameters(; data["waves"]...)
#     materials = Dict(name => Material(d["type"], d["properties"]) for (name, d) in data["materials"])
#     motions = Dict(id => Motion(d["type"], d["properties"]) for (id, d) in data["motions"])
#     seabed = SeaBedParameters(; data["seabed"]...)

#     return ParameterHandler(points, segments, lines, drag, waves, materials, motions, seabed)
# end

# function to_dict(ph::ParameterHandler)::Dict
#     return Dict(
#         "points"    => Dict(id => Dict("tag" => p.tag, "coordinates" => p.coordinates,
#                                        "motion" => p.motion, "mesh_size" => p.mesh_size) for (id,p) in ph.points),
#         "segments"  => Dict(id => Dict("tag"=>s.tag, "start"=>s.start_point,
#                                        "stop"=>s.stop_point, "length"=>s.length,
#                                        "material"=>s.material, "density"=>s.density, "area"=>s.area) for (id,s) in ph.segments),
#         "lines"     => Dict(id => Dict("points"=>l.points, "segments"=>l.segments) for (id,l) in ph.lines),
#         "drag"      => Dict(name => Dict(field=>getfield(d, field) for field in fieldnames(DragProperties)) for (name,d) in ph.drag),
#         "waves"     => Dict(field=>getfield(ph.waves, field) for field in fieldnames(WaveParameters)),
#         "materials" => Dict(name => Dict("type"=>m.type, "properties"=>m.properties) for (name,m) in ph.materials),
#         "motions"   => Dict(id => Dict("type"=>m.type, "properties"=>m.properties) for (id,m) in ph.motions),
#         "seabed"    => Dict(field=>getfield(ph.seabed, field) for field in fieldnames(SeaBedParameters))
#     )
# end

# # --- Query/Modify ---
# function get_parameter(ph::ParameterHandler, category::Symbol, key)
#     return getfield(ph, category)[key]
# end

# function set_parameter!(ph::ParameterHandler, category::Symbol, key, value)
#     getfield(ph, category)[key] = value
# end

# end # module
