module TangentialDiffCalculus

using Gridap.CellData
using Gridap.TensorValues
using Gridap.Geometry: get_triangulation, num_cell_dims

export J, G, Q, ∇ₓΓdir, FΓ, P
export j, g, q, Λ, Edir, Etang

"""
    J(X::CellField)

**Jacobian operator**

This function returns the Jacobian operator of a given coordinate map field `X`. We use the transpose as 
the gradient operator is defined as `∇(X)=∇×X` in Gridap. The Jacobian operator is defined as:

```math
mathbf{J}(mathbf{X}) = ∇_r(mathbf{X}) = rac{∂X_i}{∂r_j}
```
Note that the dimensions of the Jacobian operator are `n×1`, where `n` is the dimension of the physical space.
"""
J(X::CellField) = ∇(X)'

"""
    G(J)

**Metric tensor**

This function returns the metric tensor of a given coordinate map field `X`. The metric tensor is defined as:
\$ \\mathbf{G} = \\mathbf{J}^T \\cdot \\mathbf{J} \$, with \$ \\mathbf{J} \$ the [J](@ref) operator of the map field `X`.

Note that the dimensions of the metric tensor are `1×1`.
"""
G(J) = J'⋅J

"""
    Q(J)

**Transformation matrix**

This function returns the transformation matrix \$ \\mathbf{Q} \$ of a given coordinate map field `X`. 
The transformation matrix is defined as: 

```math
\\mathbf{Q} = \\mathbf{J}\\cdot\\mathbf{G}^{-1}
```

Note that the dimensions of the transformation matrix are `n×1`, where `n` is the dimension of the physical space.
"""
Q(J) = J⋅inv(G(J))

"""
    ∇ₓΓdir(u::CellField,X::CellField)

**Tangential gradient**

This function returns the tangential gradient of a given vector field `u`. The tangential gradient is defined as:

```math
∇_{x}^{Γ,dir}(u) = ∇_r(u)⋅Q
```

Note that the dimensions of the tangential gradient are `n×n`, where `n` is the dimension of the physical space.
"""
∇ₓΓdir(u::CellField,X::CellField) = ∇(u)'⋅(Q(J(X))')

"""
    ∇ₓΓdir(∇u::TensorValue,J::TensorValue)

[∇ₓΓdir](@ref) (Tangential gradient)
"""
∇ₓΓdir(∇u::TensorValue,J::TensorValue) = ∇u'⋅(Q(J)')

"""
    FΓ(u::CellField,X::CellField)

**Line deformation gradient**

This function returns the line deformation gradient for a given displacement field `u`, and 
a map field `X`. The line deformation gradient is defined as:

```math
F_Γ = ∇_{x}^{Γ,dir}(u) + I
```

where `I` is the identity tensor. The dimensions of the line deformation gradient are `n×n`, where `n` is the dimension of the physical space.
"""
function FΓ(u::CellField,X::CellField)
    ∇u = ∇ₓΓdir(u,X)
    I = one∘(∇u)
    return ∇u + I
end

"""
    FΓ(∇xΓdir_u::TensorValue)

[∇xΓdir](@ref) (Line deformation gradient)
"""
FΓ(∇xΓdir_u::TensorValue) = ∇xΓdir_u + one(∇xΓdir_u)

"""
    P(J)

**Projection operator**

This function returns the projection operator of a given vector field `u` onto the tangent space of the map field `X`.
The projection operator is defined as:

```math
\\mathbf{P} = \\frac{ \\mathbf{J}⋅\\mathbf{J} }{ \\|\\mathbf{J}\\|}^2}}
```

Note that the dimensions of the projection operator are `n×n`, where `n` is the dimension of the physical space.
"""
P(J) = (J⋅J')/(J⊙J)

"""
    j(FΓ,J)

**Jacobian operator in the physical space**

This function returns the Jacobian determinant in the physical space for a given line deformation gradient `FΓ` and
Jacobian operator `J`. The Jacobian determinant in the physical space is defined as:

```math
j = FΓ^T⋅J
```

Note that the dimensions of the Jacobian determinant in the physical space are `n×1`.
"""
j(FΓ,J) = FΓ'⋅J

"""
    g(j)

**Metric tensor in the physical space**

This function returns the metric tensor in the physical space as:

```math
g = j^T⋅j
```

Note that the dimensions of the metric tensor in the physical space are `1×1`.
"""
g(j) = j'⋅j

"""
    q(j)

**Transformation matrix in the physical space**

This function returns the transformation matrix in the physical space as:

```math
q = j⋅g^{-1}
```

Note that the dimensions of the transformation matrix in the physical space are `n×1`, where `n` is the dimension of the physical space.
"""
q(j) = j⋅inv(g(j))

"""
    Λ(j,J)

**Line stretch**

This function returns the stretch of a given line. The stretch is defined as:

```math
Λ = \\frac{\\det(g)}{\\det(G)}
```

where `g` is the metric tensor in the physical space, and `G` is the metric tensor in the reference space.
The result is a scalar.
"""
Λ(j,J) = (det(g(j))).^(0.5)/(det(G(J))).^(0.5)

"""
    T(J)

**Tangent vector**

This function returns the unit tangent vector of a given line. The tangent vector is defined as:
```math
\\mathbf{T} = \\frac{\\mathbf{J}}{\\|\\mathbf{J}\\|}
```
"""
# Tangent vector for scalar/tensor inputs.
# Note: G(J)[1] assumes the metric tensor is 1×1 (i.e., for lines). This is valid for 1D manifolds.
# If used for higher-dimensional manifolds, this may not be correct.
T(J) = J/sqrt(G(J)[1])

# Tangent vector for CellField inputs.
# This separate implementation is needed because Gridap handles CellField operations differently.
# The result is mathematically equivalent to the scalar/tensor case.
T(J::CellField) = J / (J'⋅J).^(0.5)

"""
    Edir(FΓ::TensorValue)

**Directional Green-Lagrange strain**

This function returns the directional Green-Lagrange strain for a given line deformation gradient `FΓ`. 
The directional Green-Lagrange strain is defined as:

```math
\\mathbf{E}_{\\text{dir}} = 0.5(\\mathbf{F}_Γ^T⋅\\mathbf{F}_Γ - \\mathbf{I})
```

where `I` is the identity tensor. The dimensions of the directional Green-Lagrange strain are `n×n`, where `n` is the dimension of the physical space.
"""
Edir(FΓ::TensorValue) = 0.5*(FΓ'⋅FΓ - one(FΓ))
function Edir(FΓ::CellField) 
    dims = num_cell_dims(get_triangulation(FΓ))
    0.5*(FΓ'⋅FΓ - one(TensorValue{dims,dims,Float64}))
end

"""
    Etang(P::TensorValue,Edir::TensorValue)

**Tangential Green-Lagrange strain**

This function returns the tangential Green-Lagrange strain for a given projection operator `P` and
directional Green-Lagrange strain `Edir`. The tangential Green-Lagrange strain is defined as:

```math
\\mathbf{E}_{\\text{tang}} = \\mathbf{P}⋅\\mathbf{E}_{\\text{dir}}⋅\\mathbf{P}
```

The dimensions of the tangential Green-Lagrange strain are `n×n`, where `n` is the dimension of the physical space.
"""
Etang(P::TensorValue,Edir::TensorValue) = P⋅Edir⋅P

end