using Pkg; Pkg.activate(".")
module dnlBasic

using DrWatson
@quickactivate "Mooring.jl"

using Printf
using Plots
using LinearAlgebra
using Roots: find_zero  
using LaTeXStrings
using CSV
using DataFrames
using XLSX
using Optim
using Statistics
using PCHIPInterpolation
using LineSearches
using BlackBoxOptim


```
Functions
=========
```

struct Schapery
  
    D0::Real
    N::Real  
    Dn::Vector{Real}
    λn::Vector{Real}
    g0::Vector{Float64}  
    g1::Vector{Float64}  
    g2::Vector{Float64} 
  
    function Schapery(
      D0::Real,
      Dn::Vector{<:Real}, 
      λn::Vector{<:Real};
      g0::Vector{<:Real} = [1.0],
      g1::Vector{<:Real} = [1.0],
      g2::Vector{<:Real} = [1.0])
  
      nDn = length(Dn)
      nλn = length(λn)
  
      N = min(nDn, nλn)
      
      new(D0, N, Dn[1:N], λn[1:N], g0, g1, g2)
    end
  end
  
  # General polynomial evaluator
  function poly_eval(σ, coeffs::Vector{Float64})
      return sum( c * σ^(i-1) for (i, c) in enumerate(coeffs) )
  end
  

# State update function
function update_step(Δt, pDn, σt1, σt0, qt0)
    g0_t1 = poly_eval(σt1, pDn.g0)
    g1_t1 = poly_eval(σt1, pDn.g1)
    g2_t1 = poly_eval(σt1, pDn.g2)
    g2_t0 = poly_eval(σt0, pDn.g2)

    temp1 = exp.(-pDn.λn * Δt)
    temp2 = (1 .- temp1) ./ (pDn.λn * Δt)

    qt1 = @. temp1 * qt0 + temp2 * (g2_t1 * σt1 - g2_t0 * σt0)
    
    Φn = @. pDn.Dn * (temp1 * qt0 - temp2 * σt0 * g2_t0)
    Φ = sum(Φn) * g1_t1

    Ψn = @. pDn.Dn - pDn.Dn * temp2
    Ψ = g0_t1 * pDn.D0 + g1_t1 * g2_t1 * sum(Ψn)

    ϵve = Ψ * σt1 - Φ

    return ϵve, qt1, Φ, Ψ
end

function get_ϵve(Δt, pDn, σin, q_init=nothing)
    T = eltype(σin)
    nTime = length(σin)
    # Fixed Dn values
    fixed_Dn = [8.249886018107235e-12, 9.015375868519326e-11, 1.118327662332956e-8, 
                1.0446117170113907e-6, 2.1632143367099333e-7, 9.472081458522793e-5, 
                7.180512827811753e-6]

    # Rebuild Schapery object with fixed Dn and dynamic g0, g1, g2
    pDn = Schapery(
        D0,
        fixed_Dn,
        λn,
        g0=pDn[nlamda+1:nlamda+4],
        g1=pDn[nlamda+5:nlamda+8],
        g2=pDn[nlamda+9:nlamda+12]
    )
    
    if isnothing(q_init)
        q = zeros(pDn.N, nTime)  
    elseif size(q_init, 2) == 1  
        q = hcat(q_init, zeros(pDn.N, nTime - 1)) 
    else  
        q = q_init
    end

    ϵve = similar(σin)
    Φ = similar(σin)
    Ψ = similar(σin)
    ϵve[1] = pDn.D0 * σin[1]

    for i in 2:nTime
        qt0 = q[:, i-1]  # Previous state
        ϵve[i], q[:, i], Φ[i], Ψ[i] = update_step(Δt, pDn, σin[i], σin[i-1], qt0)
    end

    return ϵve, q
end



# Residual function
function objective_fn(Δt, pDn, σin, ϵin, q_init=nothing)
    ϵve, _ = get_ϵve(Δt, pDn, σin, q_init)
    err = ϵve[2:end] .-= ϵin[2:end] 
    return norm(err, 2)
end

nlamda = 0
# function solve_for_Dn(λn, Δt, σin, ϵin, cInit, q_init=nothing)
#     # Define bounds for the optimization
#     n_opt = length(cInit[1:nlamda])
#     lower_bound = zeros(n_opt)
#     upper_bound = fill(1e10, n_opt)
#     initial_guess = cInit[1:nlamda]

#     # Extract the fixed values from cInit
#     fixed_g0 = cInit[nlamda+1:nlamda+4]         # g0 (locked)
#     fixed_g1 = cInit[nlamda+5:nlamda+8]         # g1 (locked)
#     fixed_g2 = cInit[nlamda+9:nlamda+12]         # g2 (locked)

#     # Define the objective function
#     function partial_objective(opt_params)
#         # Combine fixed and optimized parameters into a single vector
#         updated_params = vcat(opt_params, fixed_g0, fixed_g1, fixed_g2)
#         # Evaluate the objective function with the updated parameters
#         return objective_fn(Δt, updated_params, σin, ϵin, q_init)
#     end

#     println("Testing Partial Objective: ", partial_objective(initial_guess))

#     # Run optimization using Fminbox with Nelder-Mead
#     result = optimize(
#         partial_objective,
#         lower_bound,
#         upper_bound,
#         initial_guess,
#         Fminbox(NelderMead()),
#         Optim.Options(
#             show_trace = true,
#             iterations = 200,
#             x_tol = 1e-6
#         )
#     )

#     println("Optimization Result: ", result)

#     # Combine optimized parameters with fixed values
#     optimized_params = vcat(
#         Optim.minimizer(result),  # Optimized cInit[1:9]
#         fixed_g0,
#         fixed_g1,
#         fixed_g2
#     )

#     return optimized_params
# end


# function solve_for_Dn(λn, Δt, σin, ϵin, cInit, q_init=nothing)
#     # Define bounds for the optimization, focusing on g0, g1, g2
#     initial_guess = cInit[nlamda+1:nlamda+12]  # Initial guess for g0, g1, g2 (all 12 values)
#     n_opt = length(initial_guess)
#     lower_bound = fill(-5,n_opt)
#     upper_bound = fill(1e10, n_opt)

#     fixed_params = cInit[1:nlamda]

#     function partial_objective(opt_params)
#         updated_params = vcat(fixed_params, opt_params) 
#         return objective_fn(Δt, updated_params, σin, ϵin, q_init)
#     end

#     println("Testing Partial Objective: ", partial_objective(initial_guess))

#     # Run optimization using Fminbox with Nelder-Mead
#     result = optimize(
#         partial_objective,
#         lower_bound,
#         upper_bound,
#         initial_guess,
#         Fminbox(NelderMead()),
#         Optim.Options(
#             show_trace = true,
#             iterations = 300,
#             x_tol = 1e-6
#         )
#     )

#     println("Optimization Result: ", result)

#     # Combine optimized parameters with fixed values
#     optimized_params = vcat( 
#         fixed_params,
#         Optim.minimizer(result)
#     )

#     return optimized_params
# end

function solve_for_Dn(λn, Δt, σin, ϵin, cInit, q_init=nothing)
    # Define bounds for the optimization
    n_Dn = length(λn)
    lower_bound = vcat( fill(-1e10, 4), fill(-1e10, 4), fill(-1e10, 4))
    upper_bound = vcat( fill(1e10, 4), fill(1e10, 4), fill(1e10, 4))

    # Run optimization using Particle Swarm
    result = optimize(
        cRun -> objective_fn(Δt, cRun, σin, ϵin, q_init),  # Objective function
        lower_bound,                                       # Lower bounds
        upper_bound,                                       # Upper bounds
        cInit,                                             # Initial guess
        ParticleSwarm(                                     # Particle Swarm parameters
            n_particles = 10                               # Number of particles
        ),
        Optim.Options(                                     # Optimization options
            show_trace = true,                             # Show progress
            iterations = 1000,                             # Maximum iterations
            f_abstol = 1e-4,                               # Absolute tolerance for stopping
            f_tol = 1e-5                                   # Relative tolerance for stopping
        )
    )

    @show result

    # Return the optimized parameters
    return Optim.minimizer(result)
end



# ===================
# Main Calculations
# ===================

resDir = datadir("dnl")
fileName = resDir*"/dnl_rope2_100%_SOF2_"


# Read CSV and extract time data and experimental strain
file_path = srcdir("dnlFit", "229792_Quasi dynamic_comb.xlsx")
df = DataFrame(XLSX.readtable(file_path, "Filtered"))

df = df[2:end,:]
@show Δt_exp = df[2,3]-df[1,3] #s

## Define quantities
L0 = 10.5 #m
A = pi * (0.052/2)^2  # Cross-sectional area in m²
t = df[:,4] #time
ten = df[:,6]*1e3 #N
disp = df[:,7]/1000 #m

σ_exp = ten/A
ϵ_exp = disp/L0

# Read CSV and extract time data and experimental strain
# file_path = srcdir("dnlFit", "229793_QuasiStatic_1_all data plots.xlsx")
# df = DataFrame(XLSX.readtable(file_path, "Quasi-static"))

# df = df[2:end,:]
# @show Δt_exp = df[2,3]-df[1,3] #s

# ## Define quantities
# L0 = 10.5 #m
# A = pi * (0.052/2)^2  # Cross-sectional area in m²
# t = df[:,3] #time
# ten = df[:,5]*1e3 #N
# disp = df[:,6]/1000 #m

# σ_exp = ten/A
# ϵ_exp = disp/L0

plt1 = plot()
plot!(plt1, ϵ_exp, σ_exp, 
    lw=2, dpi=330, 
    xlabel="Strain", ylabel=" Stress (Pa)",label="Stress-Strain Experiment Data")

# get D0
ind1, ind2 = 2 , 2194
scatter!(plt1, [ϵ_exp[ind1]], [σ_exp[ind1]], label="Data Points for D0 Calculation(1)")
scatter!(plt1, [ϵ_exp[ind2]], [σ_exp[ind2]], label="Data Points for D0 Calculation(2)")
savefig(plt1, fileName*"stress_strain.png")

# @show E0 = (σ_exp[ind2]-σ_exp[ind1])/(ϵ_exp[ind2]-ϵ_exp[ind1])
# @show D0 = 1/E0
@show D0 = 1.9706432209833109e-10

# Define initial sigma and epsilon values
σ0 = σ_exp[ind1]
ϵ0 = σ0 * D0

# Update indices for coefficient definition
ind1, ind2 = 2 , 232400
scatter!(plt1, [ϵ_exp[ind1]], [σ_exp[ind1]], label="Start Points for Coefs Definition")
scatter!(plt1, [ϵ_exp[ind2]], [σ_exp[ind2]], label="End Points for Coefs Definition")
savefig(plt1, fileName * "stress_strain.png")

# Preallocate arrays with zeros
ϵr = zeros(length(ind1:ind2) + 5)
σr = zeros(length(ind1:ind2) + 5)
t_aligned = zeros(length(ind1:ind2) + 5) 

# Fill the preallocated arrays
ϵr[1:5] .= 0.0                        
ϵr[6:end] .= ϵ_exp[ind1:ind2] .+ ϵ0   
σr[1:5] .= 0.0                        
σr[6:end] .= σ_exp[ind1:ind2]        

# Fill t_aligned
t_aligned[1:5] .= t[ind1] .- (5:-1:1) .* (t[ind1] - t[ind1-1])  
t_aligned[6:end] .= t[ind1:ind2]                                


```
Run
====
```
## Define a test
λn = [1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11
] 
Dn_init = fill(1e-10, length(λn))

# cInit = vcat(Dn_init,[1],[1],[1])
cInit = [0.7273427600203776, -3.459813826638395e-10, -9.284577946815447e-18, 3.6596602310782507e-26, 0.24357602261531952, 6.634895471767583e-11, 4.289139421866027e-19, -2.6919854967211057e-27, 2.5763777426026286, -6.193118599088193e-8, 5.910734043958621e-16, -1.8674066469198892e-24]

Δt = Δt_exp
println("\nInitial error is")
@show objective_fn(Δt, cInit, σr, ϵr, nothing)


## Update After Optimization
ϵlin, qlin = get_ϵve(Δt, cInit, σr, nothing)  
@show objective_fn(Δt, cInit, σr, ϵr, nothing)

# percentage_diff = abs.((ϵr .- ϵlin) ./ ϵr) .* 100

# σCalc = zeros(length(t_aligned))
# σCalc[1] = σr[1] 

# for i in 2:length(t_aligned)
#     Δt_current = t_aligned[i] - t_aligned[i-1]
#     qt0 = qlin[:, i-1] 
#     ϵErr(σi) = ϵr[i] - get_ϵve(Δt_current, cInit, [σi, σCalc[i-1]], qt0)[1][end]
#     # Use find_zero with custom tolerances
#     σCalc[i] = find_zero(ϵErr, 0)
# end

# percentage_diff_stress = abs.((σr .- σCalc) ./ σr) .* 100

# Actual Fit

# cOpt = solve_for_Dn(λn, Δt,
#     σr, ϵr, cInit, nothing)
# @show cOpt

# ϵve = get_ϵve(Δt, cOpt, σr)
# @show objective_fn(Δt, cOpt, σr, ϵr, nothing)

# println("\nInitial error was")
# @show objective_fn(Δt, cInit, σr, ϵr, nothing)

# println("D0 used is", D0)


# -------------------- End --------------------

# ===================
# Plotting Section
# ===================

# Define the plotting paths
plt5 = plot()
plot!(plt5, t_aligned, σr, lw=2, dpi=330, 
    xlabel="Time", ylabel="Stress (Pa)", label="Time - Smooth-ed Stress Data")
# plot!(plt5, t_aligned, σCalc, lw=2, dpi=330, 
#     xlabel="Time", ylabel="Calculated Stress (Pa)", label="Inversion Stress from Modeled Strain")
savefig(plt5, fileName * "sigma_ts.png")

plt3 = plot()
plot!(plt3, t_aligned, ϵr, lw=2, dpi=330, 
    xlabel="Time", ylabel="Strain", label="Time - Smooth-ed Strain Data")
plot!(plt3, t_aligned, ϵlin, lw=2, dpi=330, 
    xlabel="Time", ylabel="Strain", label="Time - Modeled Strain Data")
savefig(plt3, fileName * "epsilon_ts.png")

plt2 = plot()
plot!(plt2, ϵr, σr, 
    lw=2, dpi=330, 
    xlabel="Strain", ylabel="Stress (Pa)", label="Smooth-ed Strain - Stress Data")
plot!(plt2, ϵlin, σr, 
    lw=2, dpi=330, 
    xlabel="Strain", ylabel="Stress (Pa)", label="Modeled Strain - Stress Data")
savefig(plt2, fileName * "trial.png")

# Plot the percentage difference for strain
# plt_diff = plot()
# plot!(plt_diff, t_aligned, percentage_diff, lw=2, dpi=330,
#     xlabel="Time", ylabel="Percentage Difference (%)",
#     label="Percentage Difference between Smooth-ed and Modeled Strain Data")
# savefig(plt_diff, fileName * "percentage_difference.png")

# Plot the percentage difference for stress
# plt_stress_diff = plot()
# plot!(plt_stress_diff, t_aligned, percentage_diff_stress, lw=2, dpi=330,
#     xlabel="Time", ylabel="Percentage Difference (%)",
#     label="Percentage Difference between Experimental and Calculated Stress")
# savefig(plt_stress_diff, fileName * "stress_percentage_difference.png")

end