#= 

Waitinng time distribution made for function calling.

=# 

# Importing Libraries 
using LinearAlgebra
using PyPlot
using Random
include("photon_counting.jl")
pygui(:qt)

# Methods
function waiting_3D_distribution(time_steps, end_time, num_of_simulations)

    Γ_list = LinRange(0.5, 5, 50)
    waiting_time_list1 = LinRange(0,end_time,1000)
    waiting_time_mat = zeros(size(Γ_list)[1],size(waiting_time_list1)[1])


    for index in 1:size(Γ_list)[1]

        Γ = Γ_list[index]
        waiting_time_list2 = return_waiting_distribution(time_steps, end_time, num_of_simulations, Γ)

        waiting_time_mat[index,:] += waiting_time_list2

        print(string(index) * "/" * string(size(Γ_list)[1]) * " sets completed", "\r")

    end

    return Γ_list, waiting_time_list1, waiting_time_mat


end

function plot_3D_waiting_distribution(Γ_list, waiting_time_list1, waiting_time_mat)

    # Remove some bits to get rid of artefacts
    waiting_time_list1 = waiting_time_list1[5:end-5]
    waiting_time_mat = waiting_time_mat[:, 5:end-5]

    surf(Γ_list, waiting_time_list1, waiting_time_mat', cmap=ColorMap("plasma"), alpha=0.8, linewidth=0.25)
    xlabel("Γ")
    ylabel("γt") 
    zlabel("w(γt)")
    plt.gca().invert_xaxis()
    gca()[:view_init](30,50)
    plt.savefig("Figures/photon_counting_3D_view12.png", dpi=600)
    gca()[:view_init](15,0)
    plt.savefig("Figures/photon_counting_3D_view22.png", dpi=600)

end

time_steps = 10000
end_time = 9
num_of_simulations = 10000

@time Γ_list, waiting_time_list1, waiting_time_mat = waiting_3D_distribution(time_steps, end_time, num_of_simulations)
plot_3D_waiting_distribution(Γ_list, waiting_time_list1, waiting_time_mat)