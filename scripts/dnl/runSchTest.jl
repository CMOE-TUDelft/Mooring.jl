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

Dn = [
    23.6e-12, 5.66e-12, 14.84e-12, 
    18.88e-12, 28.58e-12, 40.05e-12,
    60.42e-12, 79.64e-12, 162.17e-12]
λn = [1/10^(i-1) for (i,d) in enumerate(Dn)]


S = StressNLVE.Schapery(true,
    D0 = 270e-12,
    Dn = Dn,
    λn = λn,
    g0 = [1.055, -0.0007e-6, -0.0001e-12],
    g1 = [4.8, -0.45e-6, 0.017e-12, -0.0002e-18],
    g2 = [0.0, -0.0179e-6, 0.008e-12, -0.0003e-18, 3.82e-30]
  )

@show S


## Time-parameters
Δt = 1e-1
tEnd = 1500
t0 = 0
t = [t0:Δt:tEnd;]



## Strain step
function getStrain_Up(t, ϵMax)
  ϵ = t*0 .+ ϵMax
  ϵ[1] = 0

  return ϵ
end
function getStrain_UpDown(t, ϵMax)
  ϵ = ifelse.( t .< 1800, t*0 .+ ϵMax, 0.0)
  ϵ[1] = 0

  return ϵ
end
# σ = getStress_creep(t, 15e6)
ϵve = getStrain_UpDown(t, 0.010)



## Stress calc: Newton Raphson
q = zeros(S.N, length(t))
σ = zeros(length(t))
for i in firstindex(t)+1:lastindex(t)  
  
  local qt1, qt0, ϵt0, ϵt1, σt0, σt1

  qt0 = q[:,i-1]
  ϵt0 = ϵve[i-1]
  σt0 = σ[i-1]
  ϵt1 = ϵve[i]

  err(σi) = StressNLVE.Residual_σPredicted( 
    S, 
    ϵt0, Δt, qt0, σt0,
    ϵt1, Δt, σi )
  
  σt1 = find_zero(err, ϵt1/S.D0)

  
  σ[i] = σt1

  q[:,i] = StressNLVE.retqnt1.(Ref(S), λn, Ref(Δt), qt0, Ref(σt0), Ref(σt1))  

end
σNR = σ


## Stress calc: Predictor-corrector
q = zeros(S.N, length(t))
σ = zeros(length(t))
for i in firstindex(t)+1:lastindex(t)  
  
  local qt1, qt0, ϵt0, ϵt1, σt0, σt1, σtk

  qt0 = q[:,i-1]
  ϵt0 = ϵve[i-1]
  σt0 = σ[i-1]
  ϵt1 = ϵve[i]

  σtk = σt0
  for i = 1:5
    σtk, err = StressNLVE.σPredicted( 
      S, 
      ϵt0, Δt, qt0, σt0,
      ϵt1, Δt, σtk )
    # @show σtk, err
  end
  σt1 = σtk

  σ[i] = σt1

  q[:,i] = StressNLVE.retqnt1.(Ref(S), λn, Ref(Δt), qt0, Ref(σt0), Ref(σt1))  

end
σPC = σ



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
  plot!( plt, 
    xlabel = "t (s)")
end

plt1 = plot()
plt = plt1
plot!( plt, t, ϵve, 
  label = nothing, linewidth=3)
# yticks_vals = 0:0.001:0.006
# yticks_labels = [@sprintf("%.1e", y) for y in yticks_vals] 
# yticks_labels = [1e3*y for y in yticks_vals] 
# xticks_vals = 0:600:3600
# xticks_labels = xticks_vals
# plot!(plt, 
#   ylabel = L"Viscoelastic Strain $\epsilon_{ve} \times 10^{-3}$",
#   yticks = (yticks_vals,yticks_labels),
#   xticks = (xticks_vals,xticks_labels),
#   ylims=[0.0, 0.006], xlims=[-150,3600])
plot!(plt, title = "Strain Profile")
setGridlines(plt)
savefig(plt, pltName*"plt1.png")


plt3 = plot()
plt = plt3
plot!( plt, t, σNR, 
  linewidth=3, 
  label = "Newton-Raphson")
plot!( plt, t, σPC, 
  linewidth=3,
  label="Predictor-corrector")
# yticks_vals = 0:5e6:20e6
# yticks_labels = [@sprintf("%.1e", y) for y in yticks_vals] 
# yticks_labels = [1e-6*y for y in yticks_vals] 
# xticks_vals = 0:600:3600
# xticks_labels = xticks_vals
# plot!(plt, 
#   ylabel = L"Stress $\sigma$ (MPa)",
#   yticks = (yticks_vals,yticks_labels),
#   xticks = (xticks_vals,xticks_labels),
#   ylims=[0.0, 20e6], xlims=[-150, 3600])
plot!(plt, title = "Stress Profile")
setGridlines(plt)
savefig(plt, pltName*"plt3.png")



plt4 = plot()
plt = plt4
plot!( plt, ϵve, σ, 
  label = nothing, linewidth=3)
plot!(plt, title = "Stiffness Profile")
setGridlines(plt)
savefig(plt, pltName*"plt4.png")

# @show df = DataFrame(hcat(t,σCalc, ϵve), :auto)
# CSV.write("dataOut.csv",df)
