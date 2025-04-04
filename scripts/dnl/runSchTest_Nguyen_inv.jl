module testSchapery

using DrWatson
@quickactivate "Mooring.jl"

using Printf
using Plots
using Plots.PlotMeasures
using LinearAlgebra
using Roots: find_zero
using LaTeXStrings
using CSV
using DataFrames
using Mooring.StressNLVE



pltName = datadir("sims_202501","run")*"/dnl_"


## StressNLVE
# ---------------------Start---------------------   

# Nguyen (2022) Table 1
# PMMA
sch = StressNLVE.Schapery(true,
  D0 = 270.9e-12,
  Dn = [  
    23.6e-12, 5.66e-12, 14.84e-12, 
    18.88e-12, 28.58e-12, 40.05e-12,
    60.42e-12, 79.64e-12, 162.17e-12],
  λn = [  
    1.0,    1e-1,   1e-2,     
    1e-3,   1e-4,   1e-5,     
    1e-6,   1e-7,   1e-8],
  g0 = (20e6, [1.055, -7e-10, -1e-16]),
  g1 = (20e6, [4.8, -0.45e-6, 0.017e-12, -0.0002e-18]),
  g2 = (20e6, [0.0, -0.0179e-6, 0.008e-12, -0.0003e-18, 3.82e-30])
)
@show sch.g2


# # Nguyen (2022) Table 3
# sch = StressNLVE.Schapery(true,
#   D0 = 2.205e-10,
#   Dn = [  
#     2.23e-10, 2.27e-10, 1.95e-10,
#     3.50e-10, 5.50e-10, 5.50e-10],
#   λn = [  
#     1.0,    1e-1,   1e-2,     
#     1e-3,   1e-4,   1e-5],
#   g0 = [1.0],
#   g1 = [0.9875, 0.01812e-6, -3.039e-16],
#   # g2 = (0e10, [0.997, 1.493e-6, -0.0167e-12, 7e-22])
#   # g2 = [0.997, 0.0167e-6, 0.0007e-12] #NguyenCorrected
#   g2 = [1.3544, -1.8551e-8, 2.3966e-15] #Kevin
# )
# ----------------------End----------------------  

@show sch
@show sch.g0
@show sch.g1
@show sch.g2


## Time-parameters
Δt = 1e-1
tEnd = 3600
t0 = 0
t = [t0:Δt:tEnd;]


## Stress step
function getSignal_UpDown(t, ϵ1, ϵ2=0.0)
  ϵ = ifelse.( t .< 1800, ϵ1, ϵ2)
  ϵ[1] = 0

  return ϵ
end
σ = getSignal_UpDown(t, 15e6)


# ## Time-parameters
# Δt = 1e-1
# tEnd = 20
# t0 = 0
# t = [t0:Δt:tEnd;]


# ## Stress step
# function getSignal_SawTooth(t, ϵRate)  
#   tEnd = t[end]
#   tHalf = tEnd/2.0
#   ϵ = ifelse.( t .< tHalf, t.*ϵRate, (tEnd .- t).*ϵRate)  

#   return ϵ
# end
# σ = getSignal_SawTooth(t, 1e6)  


ϵve = zeros(length(t))
q = zeros(sch.N, length(t))
for i in firstindex(t)+1:lastindex(t)  
  
  local qt1, qt0, ϵt0, ϵt1, σt0, σt1

  qt0 = q[:,i-1]
  ϵt0 = ϵve[i-1]
  σt0 = σ[i-1]
  σt1 = σ[i]  

  ϵt1, err = StressNLVE.retϵve(sch, ϵt0, Δt, qt0, σt0, Δt, σt1)

  ϵve[i] = ϵt1  

  q[:,i] = StressNLVE.retqnt1.(Ref(sch), sch.λn, Ref(Δt), qt0, Ref(σt0), Ref(σt1))  

end


## Inverse function
σCalc = zeros(length(t))
q = zeros(sch.N, length(t))
for i in firstindex(t)+1:lastindex(t)  
  
  local qt1, qt0, ϵt0, ϵt1, σt0, σt1

  qt0 = q[:,i-1]
  ϵt0 = ϵve[i-1]
  ϵt1 = ϵve[i]
  σt0 = σCalc[i-1]  

  σt1, err1 = StressNLVE.get_stressNLVE(
    sch, Δt,
    ϵt0, qt0, σt0,
    ϵt1, σt0 )

  σCalc[i] = σt1

  q[:,i] = StressNLVE.retqnt1.(Ref(sch), sch.λn, Ref(Δt), qt0, Ref(σt0), Ref(σt1))  
  
end



## Plotting
function setGridlines(plt)  
  plot!(plt, grid=:true, gridcolor=:black, 
    gridalpha=0.5, gridlinestyle=:dot,
    gridlinewidth = 1)    
  plot!( plt, 
    titlefontsize=15,
    tickfontsize=13, 
    labelfontsize=15,
    legendfontsize = 13)    
  plot!(plt, size = (600, 400), dpi = 330,
    left_margin=2mm, right_margin=5mm,
    bottom_margin=2mm, top_margin=2mm)  
end

plt1 = plot()
plt = plt1
plot!( plt, t, ϵve, 
  label = nothing, linewidth=3, minorgrid=true)
plot!(plt, title = "Strain Profile")
plot!(plt, 
  xlabel = "t (s)",
  # xlims=[0,3600],
  # yticks = (0:0.001:0.015, 0:0.001:0.015)
  )
setGridlines(plt)
savefig(plt, pltName*"plt1.png")


plt2 = plot()
plt = plt2
# plot!( plt, t, σ./1e6, 
#   label = nothing, linewidth=3, minorgrid=true)
plot!( plt, t, σCalc./1e6, 
  label = nothing, linewidth=3, minorgrid=true)
plot!(plt, title = "Stress Profile MPa")
plot!(plt, 
  xlabel = "t (s)",
  # xlims=[0,1800]
  )
setGridlines(plt)
savefig(plt, pltName*"plt2.png")


plt3 = plot()
plt = plt3
plot!( plt, ϵve, σ./1e6, 
  label = nothing, linewidth=3, minorgrid=true)
plot!(plt, title = "Stress Strain Profile")
plot!(plt, 
  xlabel = L"$ \epsilon $",
  xticks = (0:0.002:0.012, 0:0.002:0.012),
  )
setGridlines(plt)
savefig(plt, pltName*"plt3.png")

CSV.write(
  pltName*"pltdata.csv", 
  DataFrame(t=t, ϵve=ϵve, σ=σ, σCalc=σCalc)
)

end 
  