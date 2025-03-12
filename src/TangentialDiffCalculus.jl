module TangentialDiffCalculus

using Gridap.CellData
using Gridap.TensorValues

export J, G, Q, ∇ₓΓdir, FΓ, P
export j, g, q, Λ

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
∇ₓΓdir(u::CellField,X::CellField) = ∇(u)'⋅(Q(J(X))')
∇ₓΓdir(∇u::TensorValue,J::TensorValue) = ∇u'⋅(Q(J)')

"""
FΓ (line deformation gradient)

This function returns the line deformation gradient for a given displacement field `u`, and 
a map field `X`. The line deformation gradient is defined as:

```math
F_Γ = ∇ₓΓdir(u) + I
```

where `I` is the identity tensor. The dimensions of the line deformation gradient are `n×n`, where `n` is the dimension of the physical space.
"""
function FΓ(u::CellField,X::CellField)
    ∇u = ∇ₓΓdir(u,X)
    I = one∘(∇u)
    return ∇u + I
end
FΓ(∇xΓdir_u::TensorValue) = ∇xΓdir_u + one(∇xΓdir_u)

"""
P (Projection operator)

This function returns the projection operator of a given vector field `u` onto the tangent space of the map field `X`.
The projection operator is defined as:

```math
\\mathbf{P} = \\frac{ \\mathbf{J}⋅\\mathbf{J} }{ \\|\\mathbf{J}\\|}^2}}
```

Note that the dimensions of the projection operator are `n×n`, where `n` is the dimension of the physical space.
"""
P(J) = (J⋅J')/(J⊙J)

j(FΓ,J) = FΓ'⋅J
g(j) = j'⋅j
q(j) = j⋅inv(g(j))

Λ(j,J) = (det(g(j))).^(0.5)/(det(G(J))).^(0.5)

end