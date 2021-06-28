#= 

Waitinng time distribution made for function calling.

=# 

# Importing Libraries 
using LinearAlgebra
using Plots 
using Random
using CSV, DataFrames
include("waiting_time.jl")

# Methods
function waiting_3D_distribution(time_steps, end_time, num_of_simulations)

    Γ_list = LinRange(0.1, 10, 15)
    waiting_time_list1 = LinRange(0,end_time,time_steps)
    waiting_time_mat = zeros(size(Γ_list)[1],size(waiting_time_list1)[1])


    for index in 1:size(Γ_list)[1]

        Γ = Γ_list[index]
        waiting_time_list2 = return_waiting_distribution(time_steps, end_time, num_of_simulations, Γ)

        waiting_time_mat[index,:] += waiting_time_list2

    end

    return Γ_list, waiting_time_list1, waiting_time_mat


end

function plot_3D_waiting_distribution(Γ_list, waiting_time_list1, waiting_time_mat)

    for row in size(waiting_time_mat)[1]

        plot(waiting_time_list1, Γ_list, waiting_time_mat[row, :])

    end


end

time_steps = 1000
end_time = 9
num_of_simulations = 1000

Γ_list, waiting_time_list1, waiting_time_mat = waiting_3D_distribution(time_steps, end_time, num_of_simulations)

CSV.write("data.csv", DataFrame(waiting_time_mat))