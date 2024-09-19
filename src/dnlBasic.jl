module dnlBasic

using Printf
using Plots
using Plots.PlotMeasures
using LinearAlgebra
using Roots: find_zero
using LaTeXStrings


struct pronySeries
  N::Real
  Dn::Vector{Real}
  λn::Vector{Real}

  function pronySeries(Dn::Vector{<:Real}, 
    λn::Vector{<:Real})

    nDn = length(Dn)
    nλn = length(λn)

    N = min(nDn, nλn)
    
    new(N, Dn[1:N], λn[1:N])
  end
end



## Visco-elastic properties
a0 = 1
g0 = 1
g1 = 1
g2 = 1

D0 = 270.9e-6 / 1e6 #Pa-1
pDn = pronySeries(
  [ 23.6358e-6/1e6, 5.6602e-6/1e6, 
    14.8405e-6/1e6, 18.8848e-6/1e6,
    28.5848e-6/1e6, 40.0569e-6/1e6,
    60.4235e-6/1e6, 79.6477e-6/1e6,
    162.1790e-6/1e6], 
  [ 1e0, 1e-1,
    1e-2, 1e-3,
    1e-4, 1e-5,
    1e-6, 1e-7,
    1e-8])



## Time-parameters
Δt = 1e-1
tEnd = 3600
t0 = 0
t = [t0:Δt:tEnd;]



## Stress
function getStress_creep(t, σMax)
  σ = t*0 .+ σMax
  σ[1] = 0

  return σ
end
function getStress_creepRelax(t, σMax)
  σ = ifelse.( t.>1800, 0.0, t*0 .+ σMax)
  σ[1] = 0

  return σ
end
# σ = getStress_creep(t, 15e6)
σ = getStress_creepRelax(t, 15e6)



## Strain
ϵve = σ*0.0



## State
function update_step(Δt, pDn, σt1, σt0, qt0)
  
  t1 = exp.(-pDn.λn * Δt)
  t2 = (1 .- t1) ./ (pDn.λn * Δt)

  qt1 = t1 .* qt0 .+ t2.*(σt1 - σt0)

  Φn = pDn.Dn .* ( t1 .* qt0 .- t2 .* σt0 )
  Φ = sum(Φn)

  Ψn = pDn.Dn .- pDn.Dn .* t2
  Ψ = D0 + sum(Ψn)

  ϵve = Ψ * σt1 - Φ

  return ϵve, qt1, Φ, Ψ
end



## Time stepping
q = zeros(pDn.N, length(t))
Φ = zeros(length(t))
Ψ = zeros(length(t))
ϵve = zeros(length(t))

for i in firstindex(t)+1:lastindex(t)  
  
  local qt1, qt0

  qt0 = q[:,i-1]
  ϵve[i], q[:,i], Φ[i], Ψ[i] = update_step(Δt, pDn, σ[i], σ[i-1], qt0)  
  
  # @printf("%i, \t %.4f", i, t[i])
  # [@printf(", \t %.4e", i) for i in q[:,i]]
  # @printf(", \t %.6e", Φ[i]) 
  # @printf(", \t %.6e", Ψ[i]) 
  # @printf(", \t %.6e", ϵve[i]) 
  # @printf("\n")  

end


## Stress calc
σCalc = zeros(length(t))
for i in firstindex(t)+1:lastindex(t)  
  
  local qt1, qt0

  qt0 = q[:,i-1]
  ϵErr(σi) = ϵve[i] - update_step(Δt, pDn, σi, σ[i-1], qt0)[1]

  σi = find_zero(ϵErr, 0)
  σCalc[i] = σi
  
  # @printf("%i, \t %.4f", i, t[i])
  # [@printf(", \t %.4e", i) for i in q[:,i]]
  # @printf(", \t %.6e", Φ[i]) 
  # @printf(", \t %.6e", Ψ[i]) 
  # @printf(", \t %.6e", σi) 
  # @printf("\n")  

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
  plot!( plt, 
    xlabel = "t (s)")
end

plt1 = plot()
plt = plt1
plot!( plt, t, ϵve, 
  label = nothing, linewidth=3)
yticks_vals = 0:0.001:0.006
# yticks_labels = [@sprintf("%.1e", y) for y in yticks_vals] 
yticks_labels = [1e3*y for y in yticks_vals] 
xticks_vals = 0:600:3600
xticks_labels = xticks_vals
plot!(plt, 
  ylabel = L"Viscoelastic Strain $\epsilon_{ve} \times 10^{-3}$",
  yticks = (yticks_vals,yticks_labels),
  xticks = (xticks_vals,xticks_labels),
  ylims=[0.0, 0.006], xlims=[-150,3600])
plot!(plt, title = "Strain Profile")
setGridlines(plt)
savefig(plt, "plt1.png")


plt2 = plot()
plt = plt2
plot!( plt, t, ϵve, 
  label = nothing, linewidth=3)
yticks_vals = 0:0.001:0.006
# yticks_labels = [@sprintf("%.1e", y) for y in yticks_vals] 
yticks_labels = [1e3*y for y in yticks_vals] 
plot!(plt, 
  ylabel = L"Viscoelastic Strain $\epsilon_{ve} \times 10^{-3}$",
  yticks = (yticks_vals,yticks_labels),
  ylims=[0.0, 0.006], xlims=[1800, 3600])
plot!(plt, title = "Relaxation Stage")
setGridlines(plt)
savefig(plt, "plt2.png")


plt3 = plot()
plt = plt3
plot!( plt, t, σCalc, 
  label = nothing, linewidth=3)
yticks_vals = 0:5e6:20e6
# yticks_labels = [@sprintf("%.1e", y) for y in yticks_vals] 
yticks_labels = [1e-6*y for y in yticks_vals] 
xticks_vals = 0:600:3600
xticks_labels = xticks_vals
plot!(plt, 
  ylabel = L"Stress $\sigma$ (MPa)",
  yticks = (yticks_vals,yticks_labels),
  xticks = (xticks_vals,xticks_labels),
  ylims=[0.0, 20e6], xlims=[-150, 3600])
plot!(plt, title = "Stress Profile")
setGridlines(plt)
savefig(plt, "plt3.png")

end