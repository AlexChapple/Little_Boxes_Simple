#= 

Waiting time distribution code. 

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


function evolve(init_condition::Array{ComplexF64}, time_list, Γ, h,)

    # coeffs_list with initial condition 
    coeffs_list::Array{Array{ComplexF64}} = [init_condition]

    # Initialise arrays for taking in σ_z, σ_+, σ_-
    σ_z_list = []
    σ_Lowering_list = []
    σ_Raising_list = []

    # photon counting lists
    photon_tracking = []

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

            # photon has been emitted so need to calculate waiting time

            c1 = 1
            c2 = 0

            push!(coeffs_list, [c1 c2 0 0 0 0]) 

            append!(σ_z_list, modulo(c2)^2 - modulo(c1)^2)
            # append!(σ_z_list, modulo(c2)^2)
            append!(σ_Lowering_list, (conj(c1)*c2))
            append!(σ_Raising_list, (c1*conj(c2)))

            append!(photon_tracking, 1)

        else

            total = sqrt(modulo(η0_new)^2 + modulo(ξ0_new)^2)

            c1 = η0_new / total
            c2 = ξ0_new / total

            push!(coeffs_list, [c1 c2 0 0 0 0]) 

            append!(σ_z_list, modulo(c2)^2 - modulo(c1)^2)
            # append!(σ_z_list, modulo(c2)^2)
            append!(σ_Lowering_list, (conj(c1)*c2))
            append!(σ_Raising_list, (c1*conj(c2)))

            append!(photon_tracking, 0)

        end


    end

    return coeffs_list, σ_z_list, σ_Lowering_list, σ_Raising_list, photon_tracking

end

function average_simulation(time_steps, end_time, num_of_simulations, Γ)

    init_condition::Array{ComplexF64} = [1 0 0 0 0 0]

    time_list = LinRange(0,end_time,time_steps)
    h = end_time/time_steps

    avg_σ_z_list = zeros(size(time_list)[1])
    avg_σ_Lowering_list = zeros(size(time_list)[1])
    avg_σ_Raising_list = zeros(size(time_list)[1])

    waiting_time_list1 = LinRange(0,end_time,1000)
    waiting_time_list2 = zeros(size(waiting_time_list1))

    for i in 1:num_of_simulations

        coeffs_list, σ_z_list, σ_Lowering_list, σ_Raising_list, photon_tracking = evolve(init_condition, time_list, Γ, h)

        # Elementwise addition of the operator lists
        avg_σ_z_list += σ_z_list
        avg_σ_Lowering_list += σ_Lowering_list
        avg_σ_Raising_list += σ_Raising_list

        # Do waiting time statistics here
        last_found_time = 0

        for q in 1:size(photon_tracking)[1]

            if photon_tracking[q] == 1 # if photon has been emitted
                waiting_time = time_list[q] - last_found_time # finds the time interval since last emission 

                for k in 1:size(waiting_time_list1)[1]

                    if waiting_time < waiting_time_list1[2]

                        waiting_time_list2[1] += 1

                    elseif waiting_time >= waiting_time_list1[end-1] && waiting_time <= waiting_time_list1[500]

                        waiting_time_list2[end] += 1

                    elseif waiting_time >= waiting_time_list1[k] && waiting_time <= waiting_time_list1[k+1]

                        waiting_time_list2[k] += 1

                    end
  

                end

                last_found_time = time_list[q]

        
            end

            

        end

        print(string(i)*"\r") # Prints to console where the code is at 

    end

    avg_σ_z_list /= num_of_simulations
    avg_σ_Lowering_list /= num_of_simulations
    avg_σ_Raising_list /= num_of_simulations

    # Return only real compoennts 
    avg_σ_z_list = [real(i) for i in avg_σ_z_list]
    avg_σ_Lowering_list = [real(i) for i in avg_σ_Lowering_list]
    avg_σ_Raising_list = [real(i) for i in avg_σ_Raising_list]

    # waiting_time_list = [i*h for i in waiting_time_list]
    waiting_time_list2 /= num_of_simulations

    return time_list, avg_σ_z_list, avg_σ_Lowering_list, avg_σ_Raising_list, waiting_time_list1, waiting_time_list2

end

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

function plot_results(time_list, avg_σ_z_list, avg_σ_Lowering_list, avg_σ_Raising_list, analytical_solutions_list, mode)

    if mode == "z"

        plot(time_list, avg_σ_z_list, lw=2, label="first order",dpi=600)
        plot!(time_list, analytical_solutions_list, label="analytic solution")
        xlabel!("\$\\gamma t\$")
        ylabel!("\$\\sigma_{z}\$")
        savefig("figure1.png")

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

function plot_waiting_time(waiting_time_list1, waiting_time_list2)

    # Needs to change dictionary to array
    plot(waiting_time_list1, waiting_time_list2,dpi=600,legend=false)
    xlabel!("γt")
    ylabel!("w(τ)")
    title!("waiting time distribution")

    savefig("waiting_max_3.5.png")
    
end

# RUN CODE

time_steps = 10000
end_time = 9
num_of_simulations = 10000
mode = "z"

Γ = 3.5

@time time_list, avg_σ_z_list, avg_σ_Lowering_list, avg_σ_Raising_list, waiting_time_list1, waiting_time_list2 = average_simulation(time_steps, end_time, num_of_simulations, Γ)
plot_waiting_time(waiting_time_list1, waiting_time_list2)



