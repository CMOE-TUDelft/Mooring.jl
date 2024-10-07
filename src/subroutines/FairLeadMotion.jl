module FairLeadMotion

using Revise
using Interpolations
using WaveSpec.Constants
using WaveSpec.WaveTimeSeries



"""
Custom Structs
=============

"""
# ---------------------Start---------------------
abstract type MotionType end
struct Regular <:MotionType end
struct WithWave <:MotionType end
# ----------------------End----------------------



"""
Functions
=============

"""
# ---------------------Start---------------------
function posItp(sp::SpecStruct, 
  inputRamp, t, x, z)  
  
  η, px, py = waveAiry1D_pPos(sp, t, x, z)      

  tRamp = timeRamp(t, inputRamp)

  itp = Interpolations.interpolate(
    η.*tRamp, 
    BSpline(Cubic(Line(OnGrid()))) )
  sitp_η = scale(itp, t)

  itp = Interpolations.interpolate(
    px.*tRamp, 
    BSpline(Cubic(Line(OnGrid()))) )
  sitp_px = scale(itp, t)

  itp = Interpolations.interpolate(
    py.*tRamp, 
    BSpline(Cubic(Line(OnGrid()))) )
  sitp_py = scale(itp, t)

  return sitp_η, sitp_px, sitp_py
end
# ----------------------End----------------------

end