module dnlBasic

using Printf
using Plots
using LinearAlgebra
using Roots: find_zero


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
Δt = 1e-2
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
  plot!( plt, dpi = 330,
    titlefontsize=15,
    tickfontsize=13, 
    labelfontsize=15,
    legendfontsize = 13)    
  plot!( plt, 
    xlabel = "t (s)")
end

plt1 = plot()
plt = plt1
plot!( plt, t, ϵve, 
  label = nothing, linewidth=3)
setGridlines(plt)
yticks_vals = 0:0.005:0.01
plot!(plt, 
  yticks = (yticks_vals,yticks_vals),
  ylims=[0.0, 0.01])
plot!(plt, xlims=[0,1800])
savefig(plt, "plt1.png")


plt2 = plot()
plt = plt2
plot!( plt, t, ϵve, 
  label = nothing, linewidth=3)
setGridlines(plt)
yticks_vals = 0:0.0005:0.001
plot!(plt, 
  yticks = (yticks_vals,yticks_vals),
  ylims=[0.0, 0.001])
plot!(plt, xlims=[1800,3600])
savefig(plt, "plt2.png")


plt3 = plot()
plt = plt3
plot!( plt, t, σCalc, 
  label = nothing, linewidth=3)
setGridlines(plt)
yticks_vals = 0:5e6:20e6
plot!(plt, 
  yticks = (yticks_vals,yticks_vals),
  ylims=[0.0, 20e6])
plot!(plt, xlims=[0,3600])
savefig(plt, "plt3.png")

end