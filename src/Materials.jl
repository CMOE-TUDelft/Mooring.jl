module Materials

using Parameters
using Gridap.TensorValues
using Gridap.CellData
import Mooring.TangentialDiffCalculus as TDC

"""
Material Struct

This abstract type is used to define the type of material. 
Possible implemented options are:
- `LinearElastic`: Elastic material
- `Scharpery`: Nonlinear viscoelastic material following Scharpery model
- `CustomMaterial`: Custom material
"""
abstract type Material end

"""
LinearElastic Struct

This struct is used to define the properties of the linear elastic material. This includes:
    - E::Real: Young's modulus [Pa]
    - μ::Real: second Lame constant [Pa]
"""
@with_kw struct LinearElastic <: Material
  E::Real = 1.0
  μ::Real = 0.5*E # second Lame constant
end

"""
Scharpery Struct

This struct is used to define the properties of the Scharpery material. This includes:
    - D0::Real: Linear elastic compliance [Pa^-1]
    - N::Int: number of relaxation times [-]
    - Dn::Vector{Real}: compliance for each relaxation time [Pa^-1] 
    - λn::Vector{Real}: Inverse of the relaxation times [s^-1]
    - g0::Function: nonlinear coefficient for the instantaneous compliance
    - g1::Function: nonlinear coefficient for the transient compliance
    - g2::Function: nonlinear coefficient for the stress rate-dependent compliance
"""
@with_kw struct Scharpery <: Material
  D0::Real = 1.0
  N::Int = 1
  Dn::Vector{Real} = [1.0]
  λn::Vector{Real} = [1.0]
  g0::Function = (σ) -> 1.0
  g1::Function = (σ) -> 1.0
  g2::Function = (σ) -> 1.0
end

function (S::Material,args...)
  error("Second Piola stress (S) not implemented for $(typeof(S)).")
end

"""
  S(material::LinearElastic,Etang)  (Second Piola stress)

This function returns the second Piola stress for a given material.
The second Piola stress is defined in terms of the derivative of the strain energy density 
with respect to the tangential Lagrangian Green Strain: \$ \\mathbf{S} = \\frac{∂W}{∂\\mathbf{Etang}} \$.

For linear elastic material this is: \$ \\mathbf{S} = 2\\mu \\mathbf{Etang} \$.
"""
S(material::LinearElastic,Etang) = 2*material.μ*Etang

"""
  S(material::Material,Xh::CellField,uh::CellField)  (Second Piola stress)

This function returns the second Piola stress for a given material. This is a generic
implementation that takes the material, the deformation map (`Xh`) and the displacement field (`uh`) as arguments.
"""
function S(material::Material,Xh::CellField,uh::CellField)
  Edir = TDC.Edir∘(TDC.FΓ(uh,Xh))
  P = TDC.P(TDC.J(Xh))
  Etang(x) = (TDC.Etang∘(P,Edir))(x)
  return x -> S(material,Etang(x))
end

"""
    K(FΓ::TensorValue,S::TensorValue)

K(FΓ::TensorValue,S::TensorValue)  (First Piola stress)

This function returns the first Piola stress for a given material. 
The first Piola stress is defined in terms of the deformation along the line (\$ FΓ \$) and 
the second Piola stress (\$ \\mathbf{S} \$) as:

```math
\\mathbf{K} = \\mathbf{FΓ}⋅\\mathbf{S}
```
"""
K(FΓ::TensorValue,S::TensorValue) = FΓ ⋅ S

"""
    σ(Λ::Real,FΓ::TensorValue,K::TensorValue)

**Cauchy stress**

This function returns the Cauchy stress for a given material.
The Cauchy stress is defined as: 

```math
\\mathbf{σ} = \\frac{1}{Λ}\\mathbf{K}⋅\\mathbf{F}_Γ^T
```

where \$ Λ \$ is the stretch of the line, \$ \\mathbf{F}_Γ \$ is the deformation gradient along the line,
and \$ \\mathbf{K} \$ is the first Piola stress.
"""
σ(Λ::Real,FΓ::TensorValue,K::TensorValue) = (1/Λ)*(K ⋅ FΓ')


end