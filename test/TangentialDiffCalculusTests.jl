module TangentialDiffCalculusTests

import Mooring.TangentialDiffCalculus as TDC
using Gridap

# Construct reference map CellField
model = CartesianDiscreteModel((0,1),(1,))
reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},1)
V = TestFESpace(model,reffe)
X(r) = VectorValue(2r[1],0.0)
Xₕ = interpolate(X,V)

end