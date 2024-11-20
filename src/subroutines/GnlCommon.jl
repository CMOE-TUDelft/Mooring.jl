module GnlCommon

using Revise
using Gridap
using Gridap.Algebra
using Gridap.ODEs
using Gridap.Arrays: testitem, return_cache
using Roots: find_zero
using DataFrames:DataFrame
using DataFrames:Matrix
using WriteVTK
using TickTock
using Parameters
using LineSearches: BackTracking
using LineSearches: Static
using LinearAlgebra
using CSV
using Interpolations
using Printf
using WaveSpec.Constants
using WaveSpec.Jonswap
using WaveSpec.WaveTimeSeries
using WaveSpec.Currents
using Plots



export printTerAndFile, showTerAndFile
export getInputSpec, getCurrentField, getWaveVel
export setInitXZ, getParabola
export assemble_cache



"""
Aux functions
=============

"""
# ---------------------Start---------------------

function printTerAndFile(str::String, 
    val::Union{AbstractArray,Tuple,Real}, outFile::IOStream)
  
  println(str, val)
  println(outFile, str, val)
  # @printf("%s %15.6f\n", str, val)
  # @printf(outFile,"%s %15.6f\n", str, val)
end

function printTerAndFile(str::String, outFile::IOStream)
  println(str)
  println(outFile, str)
end


# function showTerAndFile(varName::Symbol, var::Any, outFile::IOStream)
#   @show var
#   println(outFile, varName)
#   show(outFile, var)
#   println(outFile)
# end
# showTer(varName::Symbol, var::Any) = 
  #   showTerAndFile(varName,var,outFile0)
function showTerAndFile(var::Any, outFile::IOStream)
  @show var
  show(outFile, var)
  println(outFile)
end


function setInitXZ(initCSV)
  
  initXZ = CSV.read(initCSV, DataFrame, header=false)
  initXZ = Matrix(initXZ)

  dx = initXZ[2:end,:].-initXZ[1:end-1,:]
  dx = [0 0; dx]

  ds = zeros(size(dx,1))
  r = zeros(size(dx,1))

  for i in axes(dx,1)
    ds[i] = norm(dx[i,:])
  end

  for i in 2:length(ds)
    r[i] = r[i-1] + ds[i]
  end

  r = r / r[end]

  interpX = linear_interpolation(r, initXZ[:,1])
  interpZ = linear_interpolation(r, initXZ[:,2])
  
  # X(r) = VectorValue(interpX(r), interpZ(r))
  # Xh = interpolate_everywhere(X, Ψu)

  # return Xh

  return interpX, interpZ

end


function getParabola(xend,zend,L)

  dx = xend /100
  x = 0:dx:xend

  b(a) = zend/xend - a*xend  
  yy(a) = 2*a.*x .+ b(a)
  y(a) = sqrt.(1.0 .+ yy(a).*yy(a))

  lineLen(a) = gaussQuad1D(y(a), dx)
  errLen(a) = lineLen(a) - L
  
  a = find_zero(errLen, 0.1)  

  pA = a
  pB = b(a)  

  X(r) = VectorValue( 
    r[1]/L*xend, 
    pA * (r[1]/L*xend)^2 + pB*(r[1]/L*xend) )

  
  printMsg = @sprintf("[VAL] (pA, pB) = %15.6e, %15.6e", pA, pB)
  
  return X, printMsg
end


function assemble_cache(xNew, save_f_cache2)

  cell_f = get_array(xNew)
  cell_f_cache = array_cache(cell_f)    
  cache2 = cell_f_cache, save_f_cache2, cell_f, xNew

  return cache2

end
# ----------------------End----------------------



"""
Wave and Current functions
=============

"""
# ---------------------Start---------------------
function getInputSpec(params)

  @unpack Hs, Tp, h0, nω, seed, ωc, enableWaveSpec = params

  if(enableWaveSpec)
    if(ωc < 0)
      ω, S, A = jonswap(Hs, Tp,
        plotflag=false, nω = nω)
    else
      ω, S, A = jonswap(Hs, Tp,
        plotflag=false, nω = nω, ωc = ωc)
    end
  
  else
    ω, S, A = jonswap(0.0, 10.0,
        plotflag=false, nω = 64)
  end

  k = dispersionRelAng.(h0, ω; msg=false)
  α = randomPhase(ω, seed = seed)

  sp = SpecStruct( h0, ω, S, A, k, α; Hs = Hs, Tp = Tp )
  return sp
end


function getCurrentField(r, Xh, curObj)
  
  X_qp = Xh(r)
  pz = X_qp ⋅ VectorValue(0.0,1.0) - curObj.h0
  pz = min(pz, 0.0)
  pz = max(pz, -curObj.h0)

  return VectorValue( curObj.itp( pz ), 0.0 )   

end


# function getWaveVelField(r, t, sp, xh)

#   x_qp = xh(r)  
#   w_u, w_w = waveAiry1D_vel(sp, t, 
#     x_qp[1], x_qp[2]-sp.h0 )

#   return VectorValue(w_u, w_w)
# end


function getWaveVel(t, sp, x)
    
  w_u, w_w = waveAiry1D_vel(sp, t, 
    x[1], x[2]-sp.h0 )

  return VectorValue(w_u, w_w)
end
# ----------------------End----------------------


end