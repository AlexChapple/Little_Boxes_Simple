#= 

Waiting time distribution code. This is for single plots only. 

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

function evolve(init_condition::Array{ComplexF64}, time_list, Γ, h)

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

    return photon_tracking

end

function average_simulation(time_steps, end_time, num_of_simulations, Γ)

    init_condition::Array{ComplexF64} = [1 0 0 0 0 0]

    time_list = LinRange(0,end_time,time_steps)
    h = end_time/time_steps

    waiting_time_list1 = LinRange(0,end_time,1000)
    waiting_time_list2 = zeros(size(waiting_time_list1))

    for i in 1:num_of_simulations

        photon_tracking = evolve(init_condition, time_list, Γ, h)

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

        if i % 10 == 0
            print(string(i) * "/" * string(num_of_simulations) * " Simulations Completed.\r")
        end

    end

    # waiting_time_list = [i*h for i in waiting_time_list]
    waiting_time_list2 /= num_of_simulations

    return time_list, waiting_time_list1, waiting_time_list2

end

function plot_waiting_time(waiting_time_list1, waiting_time_list2, Γ)

    # Needs to change dictionary to array
    plot(waiting_time_list1, waiting_time_list2,dpi=600,legend=false)
    xlabel!("γt")
    ylabel!("w(τ)")
    title!("waiting time distribution")
    title = "Figures/waiting_distribution_" * string(Γ) * ".png"
    savefig(title)
    
end

# RUN CODE

time_steps = 10000
end_time = 9
num_of_simulations = 10000

Γ = 0.8

@time time_list, waiting_time_list1, waiting_time_list2 = average_simulation(time_steps, end_time, num_of_simulations, Γ)
plot_waiting_time(waiting_time_list1, waiting_time_list2, Γ)



