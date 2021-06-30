#= 

Waitinng time distribution made for function calling.

=# 

# Importing Libraries 
using LinearAlgebra
using Plots 
using Random


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

    photon_counting_time_list = LinRange(0,end_time,1000) # changed 1000 to time_steps but might cause issues
    photon_counting_list = zeros(size(photon_counting_time_list))

    for i in 1:num_of_simulations

        photon_tracking = evolve(init_condition, time_list, Γ, h)

        for q in 1:size(photon_tracking)[1]

            if photon_tracking[q] == 1 # if photon is emitted

                time = time_list[q] # Finds the time it has been found 

                if time < photon_counting_time_list[2]
                    photon_counting_list[1] += 1
                elseif time >= photon_counting_time_list[end-1] && time <= photon_counting_time_list[end]
                photon_counting_list[end] += 1
                else
                    for k in 2:(size(photon_counting_time_list)[1]-1)
                        if time >= photon_counting_time_list[k] && time < photon_counting_time_list[k+1]
                            photon_counting_list[k] += 1
                            break
                        end
                    end 
                end
            end
        end

    end

    # waiting_time_list = [i*h for i in waiting_time_list]
    photon_counting_list /= num_of_simulations

    return time_list, photon_counting_time_list, photon_counting_list

end

function return_waiting_distribution(time_steps, end_time, num_of_simulations, Γ)

    time_list, photon_counting_time_list, photon_counting_list = average_simulation(time_steps, end_time, num_of_simulations, Γ)

    return photon_counting_list

end