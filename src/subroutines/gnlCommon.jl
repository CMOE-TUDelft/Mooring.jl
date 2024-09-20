module gnlCommon

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



export printTerAndFile, getInputSpec, setInitXZ
export getParabola, assemble_cache



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


function getInputSpec(params)

  @unpack Hs, Tp, h0, nω, seed, ωc = params

  if(ωc < 0)
    ω, S, A = jonswap(Hs, Tp,
      plotflag=false, nω = nω)
  else
    ω, S, A = jonswap(Hs, Tp,
      plotflag=false, nω = nω, ωc = ωc)
  end

  k = dispersionRelAng.(h0, ω; msg=false)
  α = randomPhase(ω, seed = seed)

  sp = SpecStruct( h0, ω, S, A, k, α; Hs = Hs, Tp = Tp )
  return sp
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

  return a,b(a)
end


function assemble_cache(xNew, save_f_cache2)

  cell_f = get_array(xNew)
  cell_f_cache = array_cache(cell_f)    
  cache2 = cell_f_cache, save_f_cache2, cell_f, xNew

  return cache2

end
# ----------------------End----------------------

end