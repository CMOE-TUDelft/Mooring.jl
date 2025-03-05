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
    σt1_n = σt1#/MBL
    σt0_n = σt0#/MBL
    g0_t1 = poly_eval(σt1_n, pDn.g0)
    g1_t1 = poly_eval(σt1_n, pDn.g1)
    g2_t1 = poly_eval(σt1_n, pDn.g2)
    g2_t0 = poly_eval(σt0_n, pDn.g2)
    
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

function get_ϵve(tin, pDn, σin, q_init=nothing)
    T = eltype(σin)
    param_guess = pDn
    nTime = length(tin)
    pDn = Schapery(
    D0,
    param_guess[1:nlamda],           
    λn,
    g0=param_guess[nlamda+1:nlamda+4],       
    g1=param_guess[nlamda+5:nlamda+8],      
    g2=param_guess[nlamda+9:nlamda+12]         
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
    ϵve[1:5] .= 4*ϵ0#pDn.D0 * σin[1]
    
    for i in 6:nTime
        qt0 = q[:, i-1]  # Previous state
        Δt = tin[i]-tin[i-1]
        if Δt≈0.0; Δt=0.5 end
        ϵve[i], q[:, i], Φ[i], Ψ[i] = update_step(Δt, pDn, σin[i], σin[i-1], qt0)
    end
    
    return ϵve, q
end

function compute_strain(coeffs)
    tin = vcat(-2.0:0.5:0,df[:,4])
    σ_exp = vcat(fill(tension[1]/A,5),tension/A)

    # Compute model strains
    ϵve,qve = get_ϵve(tin, coeffs, σ_exp)

    # Reshape quantities
    ϵve = ϵve[6:end].-ϵ0
    tin = tin[6:end]
    σ_exp = σ_exp[6:end]
    
    # Strain Error
    error = abs.((ϵ_exp .- ϵve) ./ ϵ_exp) .* 100

    return tin, ϵve, σ_exp, error
end


# ===================
# Main Calculations
# ===================

resDir = datadir("dnl")
fileName = resDir*"/dnl_rope2_final_"

# Read CSV and extract time data and experimental strain
file_path = srcdir("dnlFit", "Quasi-static_rope229793.xlsx")
df = DataFrame(XLSX.readtable(file_path, "Sheet3"))
df = df[2:end,:]

## Define quantities
L0 = 10.5 #m
A = pi * (0.052/2)^2  # Cross-sectional area in m²
time = df[:,4] # time
tension = df[:,6]*1e3 #N
displacement = df[:,7]/1000 #m
strain = df[:,13]/100
@show D0 = 1.97e-10 # From quasi-static
@show ϵ0 = 1.5e-3 # Gives good approximation  
ϵ_exp = strain
#ϵ_exp .+= ϵ0
    
    

## Define parameters
λn = [100,10,1.0, 1.0e-1, 1e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8,1.0e-9,1.0e-10
]
# λn = [1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8,1.0e-9,1.0e-10
# ] 
nlamda = length(λn)
# coeffs = [1.0e-10, 5.0e-11, 1.0e-11, 5.0e-11, 1.0e-10, 1.0e-10, 2.0e-10, 9.0e-10, 1.0e-11,
# 0.50, -2.0e-10, -1.0e-17, 3.7e-26, 
# 0.24, 6.6e-11, 4.2e-19, -2.7e-27, 
# 1.2, -1.05e-8, 1.1e-16, -3.8e-25]

coeffs = [
       1.0e-10, # 1e2 
       5.0e-11, # 1e1 
       2.0e-11, # 1e0 
       8.0e-11, # 1e-1
       1.0e-10, # 1e-2
       1.0e-10, # 1e-3 
       2.0e-10, # 1e-4
       9.0e-10, # 1e-5
       1.0e-10, # 1e-6
       1.0e-10, # 1e-7
       1.0e-16, # 1e-8 
       1.0e-16, # 1e-9 
       1.0e-16,  # 1e-10
       0.50, -2.0e-10, -1.0e-17, 3.7e-26, # g0
       0.24, 6.6e-11, 4.2e-19, -2.7e-27,  # g1
       1.4, -9.8e-9, 7.4e-17, -2.0e-25]   # g2

# # 19-01-2025 (2)
# coeffs = [       1.0e-10, # 1e2 
# 5.0e-11, # 1e1 
# 2.0e-11, # 1e0 
# 8.0e-11, # 1e-1
# 1.0e-10, # 1e-2
# 1.0e-10, # 1e-3 
# 2.0e-10, # 1e-4
# 1.0e-9, # 1e-5
# 1.6e-8, # 1e-6
# 1.0e-10, # 1e-7
# 1.0446117170113907e-17, # 1e-8 
# 2.1632143367099333e-17, # 1e-9 
# 9.472081458522793e-17,  # 1e-10
# 0.7273427600203776, -3.459813826638395e-10, -9.284577946815447e-18, 3.6596602310782507e-26, 
# 0.245, -2.1e-10, 1.6e-18, -2.6919854967211057e-27,  # g1
# 1.92, -6.8e-8, 8.7e-16, -3.2e-24]   # g2

# coeffs = [
# 1.0e-10, # 1e2 
# 5.0e-11, # 1e1 
# 2.0e-11, # 1e0 
# 8.0e-11, # 1e-1
# 1.0e-10, # 1e-2
# 2.5e-10, # 1e-3 
# 3.0e-10, # 1e-4
# 8.0e-10, # 1e-5
# 6.0e-9, # 1e-6
# 1.0e-10, # 1e-7
# 0.46, 7.5e-10, -1.0e-17, 1.5e-26, # g0
# 0.24, 6.2e-10, -8.9e-18, 4.4e-26,  # g1
# 0.45, -5.0e-10, 5.0e-17, -3.0e-25]   # g2

# 1.0e-16,1.0e-16,.0e-16,1.0e-16,.0e-16,1.0e-16,1.0e-16,
# coeffs = [
#     8.249886018107235e-12, # 1e-5
#     9.015375868519326e-11, # 1e-6
#     1.118327662332956e-8,  # 1e-7 
#     1.0446117170113907e-6, # 1e-8 
#     2.1632143367099333e-7, # 1e-9 
#     9.472081458522793e-5,  # 1e-10
#     0.7273427600203776, -3.459813826638395e-10, -9.284577946815447e-18, 3.6596602310782507e-26, # g0
#     0.24357602261531952, 6.634895471767583e-11, 4.289139421866027e-19, -2.6919854967211057e-27, # g1
#     2.5763777426026286, -6.193118599088193e-8, 5.910734043958621e-16, -1.8674066469198892e-24   # g2
# ]

# Compute strains
tin, ϵve, σ_exp, error = compute_strain(coeffs)

# Postprocess
function compute_and_plot(coeffs,range=nothing)
    tin, ϵve, σ_exp, error = compute_strain(coeffs)

    plt1 = plot()
    plot!(plt1, ϵ_exp.-ϵ_exp[45], σ_exp, 
        lw=2, dpi=330, 
        xlabel="Strain", ylabel=" Stress (Pa)",label="Stress-Strain Experiment Data"
    )
    savefig(plt1, fileName * "strain_stress_experiment.png")

    plt2 = plot()
    plot!(plt2, tin, σ_exp, lw=2, dpi=330, 
    xlabel="Time", ylabel="Stress (Pa)", label="Time - Smooth-ed Stress Data")
    savefig(plt2, fileName * "sigma_ts.png")

    plt3 = plot()
    if range==nothing
        a1=1
        a2=length(tin)
    else
        a1,a2=range
    end
    # plot!(plt3, tin[a1:a2], ϵ_exp[a1:a2].-ϵ_exp[10000], lw=2, dpi=330, 
    # xlabel="Time", ylabel="Strain", label="Time - Smooth-ed Strain Data")
    # plot!(plt3, tin[a1:a2], ϵve[a1:a2].-ϵve[10000], lw=1, dpi=330, 
    # xlabel="Time", ylabel="Strain", label="Time - Modeled Strain Data")
    plot!(plt3, tin[a1:a2], ϵ_exp[a1:a2].-ϵ_exp[45], lw=0.8, dpi=330, 
    xlabel="Time", ylabel="Strain", label="Time - Smooth-ed Strain Data")
    plot!(plt3, tin[a1:a2], ϵve[a1:a2].-ϵve[45], lw=0.5, dpi=330, 
    xlabel="Time", ylabel="Strain", label="Time - Modeled Strain Data")
    savefig(plt3, fileName * "epsilon_ts.png")

    plt4 = plot()
    plot!(plt4, ϵ_exp[a1:a2].-ϵ_exp[45], σ_exp[a1:a2], 
    lw=0.8, dpi=330, 
    xlabel="Strain", ylabel="Stress (Pa)", label="Smooth-ed Strain - Stress Data")
    plot!(plt4, ϵve[a1:a2].-ϵve[45], σ_exp[a1:a2], 
    lw=0.5, dpi=330, 
    xlabel="Strain", ylabel="Stress (Pa)", label="Modeled Strain - Stress Data")
    savefig(plt4, fileName * "strain_stress_model.png")

    # Plot the percentage difference for strain
    plt_diff = plot()
    plot!(plt_diff, tin[a1:a2], error[a1:a2], lw=2, dpi=330,
    xlabel="Time", ylabel="Percentage Difference (%)",
    label="Percentage Difference between Smooth-ed and Modeled Strain Data")
    savefig(plt_diff, fileName * "error.png")

    # Write to xlsx
    XLSX.openxlsx(srcdir("dnlFit", "229792_model_data.xlsx"), mode="w") do xf
        sheet = xf[1]
        columns = [tin,ϵ_exp .* 100,ϵve .* 100]
        labels = ["time","Δl/l-exp","Δl/l-model"]
        XLSX.writetable!(sheet, columns, labels)
    end
end

compute_and_plot(coeffs)

end