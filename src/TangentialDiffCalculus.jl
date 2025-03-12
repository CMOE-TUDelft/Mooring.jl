module TangentialDiffCalculus

using Gridap.CellData

export J, G, Q, ∇ₓΓdir

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

Note that the dimensions of the metric tensor are `1×1`.
"""
G(J) = J'⋅J

"""
Q (Transformation matrix)

This function returns the transformation matrix \$ Q \$ of a given coordinate map field `X`. 
The transformation matrix is defined as: 

```math
\\mathbf{Q} = \\mathbf{J}\\cdot\\mathbf{G}^{-1}
```

Note that the dimensions of the transformation matrix are `n×1`, where `n` is the dimension of the physical space.
"""
Q(J) = J⋅inv(G(J))

"""
∇ₓΓdir (Tangential gradient)

This function returns the tangential gradient of a given vector field `u`. The tangential gradient is defined as:

```math
∇ₓΓdir(u) = ∇_r(u)⋅Q(J(X))
```

Note that the dimensions of the tangential gradient are `n×n`, where `n` is the dimension of the physical space.
"""
∇ₓΓdir(u,X) = ∇(u)'⋅(Q(J(X))')

end