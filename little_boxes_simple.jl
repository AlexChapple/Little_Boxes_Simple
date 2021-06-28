#= 

Little boxes simple with first order 
approximation of the 6 differential equations. 

=#

# Importing Libraries 
using LinearAlgebra
using Plots 
using Random

# Constants

# Methods 
function modulo(z)
    return sqrt(real(z)^2 + imag(z)^2)
end

# Evolves one simulation 
function evolve(init_condition::Array{ComplexF64}, time_list, Γ, h,)

    # coeffs_list with initial condition 
    coeffs_list::Array{Array{ComplexF64}} = [init_condition]

    # Initialise arrays for taking in σ_z, σ_+, σ_-
    σ_z_list = []
    σ_Lowering_list = []
    σ_Raising_list = []

    for t in time_list

        coeffs = coeffs_list[end] # Grabs the latest coefficients

        η0_new = coeffs[1] - (Γ/2)*coeffs[2]*h
        ξ0_new = coeffs[2]*(1 - h/2) + (Γ/2)*coeffs[1]*h
        ηLR_new = -1im * sqrt(h/2) * coeffs[2]

        total = sqrt(modulo(η0_new)^2 + modulo(ξ0_new)^2 + modulo(ηLR_new)^2)
        
        η0_new / total
        ξ0_new / total
        ηLR_new / total

        # Generate random number 
        random_num = rand()

        if random_num <= 2 * modulo(ηLR_new)^2

            c1 = 1
            c2 = 0

            push!(coeffs_list, [c1 c2 0 0 0 0]) 

            append!(σ_z_list, modulo(c2)^2 - modulo(c1)^2)
            # append!(σ_z_list, modulo(c2)^2)
            append!(σ_Lowering_list, (conj(c1)*c2))
            append!(σ_Raising_list, (c1*conj(c2)))

        else

            total = sqrt(modulo(η0_new)^2 + modulo(ξ0_new)^2)

            c1 = η0_new / total
            c2 = ξ0_new / total

            push!(coeffs_list, [c1 c2 0 0 0 0]) 

            append!(σ_z_list, modulo(c2)^2 - modulo(c1)^2)
            # append!(σ_z_list, modulo(c2)^2)
            append!(σ_Lowering_list, (conj(c1)*c2))
            append!(σ_Raising_list, (c1*conj(c2)))

        end


    end

    return coeffs_list, σ_z_list, σ_Lowering_list, σ_Raising_list

end

# Average simulations 
function average_simulation(time_steps, end_time, num_of_simulations, Γ)

    init_condition::Array{ComplexF64} = [1 0 0 0 0 0]

    time_list = LinRange(0,end_time,time_steps)
    h = end_time/time_steps

    avg_σ_z_list = zeros(size(time_list)[1])
    avg_σ_Lowering_list = zeros(size(time_list)[1])
    avg_σ_Raising_list = zeros(size(time_list)[1])

    for i in 1:num_of_simulations

        coeffs_list, σ_z_list, σ_Lowering_list, σ_Raising_list = evolve(init_condition, time_list, Γ, h)

        # Elementwise addition of the operator lists
        avg_σ_z_list += σ_z_list
        avg_σ_Lowering_list += σ_Lowering_list
        avg_σ_Raising_list += σ_Raising_list

        # Prints to console the number of simulations completed
        if i % 10 == 0
            print(string(i) * "/" * string(num_of_simulations) * " Simulations Completed.\r")
        end

    end

    avg_σ_z_list /= num_of_simulations
    avg_σ_Lowering_list /= num_of_simulations
    avg_σ_Raising_list /= num_of_simulations

    # Return only real compoennts 
    avg_σ_z_list = [real(i) for i in avg_σ_z_list]
    avg_σ_Lowering_list = [real(i) for i in avg_σ_Lowering_list]
    avg_σ_Raising_list = [real(i) for i in avg_σ_Raising_list]

    return time_list, avg_σ_z_list, avg_σ_Lowering_list, avg_σ_Raising_list

end

# Calculates analytical solution 
function analytical_solution(Γ, time_list, mode)

    solution_list = []

    Y = sqrt(2) * Γ
    δ = (1/4)*sqrt(0im + 1 - 8*Y^2)

    if mode == "z"
    
        for t in time_list

            a = -1/(1 + Y^2)
            b = Y^2 * exp(-3*t/4)
            c = cosh(δ*t)
            d = (3 / sqrt(0im + 1 - 8*Y^2)) * sinh(δ*t)

            result = a*(1 + b*(c + d))

            append!(solution_list, result)

        end

        solution_list = [real(i) for i in solution_list]

    elseif mode == "-"

        for t in time_list

            a = (1im / sqrt(2)) * (Y / (1 + Y^2)) 
            b = exp(-3t/4) 
            c = cosh(δ*t) + (3/sqrt(1 - 8*Y^2 + 0im))*sinh(δ*t)
            d = 1im*sqrt(2)*Y*b*(1/sqrt(1 - 8*Y^2 + 0im))*sinh(δ*t)

            result = a*(1 - b*c) + d
            append!(solution_list, result)

        end

        solution_list = [imag(i) for i in solution_list]

    elseif mode == "+"

        for t in time_list

            a = -(1im / sqrt(2)) * (Y / (1 + Y^2)) 
            b = exp(-3t/4) 
            c = cosh(δ*t) + (3/sqrt(1 - 8*Y^2 + 0im))*sinh(δ*t)
            d = -1im*sqrt(2)*Y*b*(1/sqrt(1 - 8*Y^2 + 0im))*sinh(δ*t)

            result = a*(1 - b*c) + d
            append!(solution_list, result)

        end

        solution_list = [imag(i) for i in solution_list]

    end

    return solution_list


end

# Plotting function 
function plot_results(time_list, avg_σ_z_list, avg_σ_Lowering_list, avg_σ_Raising_list, analytical_solutions_list, mode, Γ)

    if mode == "z"

        plot(time_list, avg_σ_z_list, lw=2, label="first order",dpi=600)
        plot!(time_list, analytical_solutions_list, label="analytic solution")
        xlabel!("\$\\gamma t\$")
        ylabel!("\$\\sigma_{z}\$")
        title = "Figures/sigma_z_" * string(Γ) * ".png"
        savefig(title)

    elseif mode == "+"

        plot(time_list, avg_σ_Raising_list, lw=2, label="first order",dpi=600)
        plot!(time_list, analytical_solutions_list, label="analytic solution")
        xlabel!("\$\\gamma t\$")
        ylabel!("\$\\sigma_{+}\$")
        savefig("figure2.png")

    elseif mode == "-"

        plot(time_list, avg_σ_Lowering_list, lw=2, label="first order",dpi=600)
        plot!(time_list, analytical_solutions_list, label="analytic solution")
        xlabel!("\$\\gamma t\$")
        ylabel!("\$\\sigma_{-}\$")
        savefig("figure3.png")

    end

end

# RUN CODE

time_steps = 10000
end_time = 20
num_of_simulations = 10000
mode = "z"

Γ = 2.3

@time time_list, avg_σ_z_list, avg_σ_Lowering_list, avg_σ_Raising_list = average_simulation(time_steps, end_time, num_of_simulations, Γ)
analytical_solutions_list = analytical_solution(Γ, time_list, mode)

plot_results(time_list, avg_σ_z_list, avg_σ_Lowering_list, avg_σ_Raising_list, analytical_solutions_list, mode, Γ)


