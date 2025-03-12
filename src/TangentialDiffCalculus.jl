module TangentialDiffCalculus

using Gridap.CellFields

"""
J (Jacobian operator)

This function returns the Jacobian operator of a given coordinate map field `X`. We use the transpose as 
the gradient operator is defined as `∇(X)=∇×X` in Gridap. The Jacobian operator is defined as:

```math
\\mathbf{J}(\\mathbf{X}) = ∇(\\mathbf{X}) = \\frac{∂X_i}{∂r_j}
```
Note that the dimensions of the Jacobian operator are `n×1`, where `n` is the dimension of the physical space.
"""
J(X::CellField) = ∇(X)'

"""
G (Metric tensor)

This function returns the metric tensor of a given coordinate map field `X`. The metric tensor is defined as:
\$ \\mathbf{G} = \\mathbf{J}^T \\cdot \\mathbf{J} \$, with \$ \\mathbf{J} \$ the [J](@ref) operator of the map field `X`.
"""
G(J) = J'⋅J

"""
Ginv (Inverse metric tensor)

This function returns the inverse of the metric tensor [G](@ref).
"""
Ginv(G) = inv(G)


end