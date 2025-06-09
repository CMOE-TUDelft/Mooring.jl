import Mooring.Drag as Drag
using Parameters

# Testing DragProperties
custom_type = Drag.Custom()
@test_throws "Invalid DragType for default constructor." default_drag_prop = Drag.DragProperties(custom_type)
custom_drag_prop = Drag.DragProperties(custom_type, 1.0, 2.0)
@test custom_drag_prop.Cfd_n == 0.5 * 1025 / 2.0

# Testing drag force
drag_force = Drag.drag_ΓX(custom_drag_prop, 1.0, VectorValue(1.0,0.0), VectorValue(1.0,0.0))
@test drag_force == VectorValue(custom_drag_prop.Cfd_t,0.0)



