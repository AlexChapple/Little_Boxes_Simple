using LinearAlgebra
using Plots

function analytical_array(end_time, time_steps, β, Ω)

    time_list = LinRange(0,end_time,time_steps)

    analytical_solution = zeros(Float64, size(time_list)[1])

    for index in 1:size(time_list)[1]
        analytical_solution[index] = real(analytical_distribution_equation(time_list[index], β, Ω))
    end  

    return time_list, analytical_solution

end

function analytical_distribution_equation(t, β, Ω)

    a = (2 * β * Ω^2) / (β^2 - Ω^2)
    b = exp(-β*t)
    c = (sinh(0.5*(sqrt(β^2 - Ω^2 + 0im))*t))^2

    return a*b*c

end

function plot_analytical_solution(time_list, analytical_solution, β, Ω)

    plot(time_list, analytical_solution, dpi=600, legend=false)
    xlabel!("γt")
    ylabel!("w(γt)")
    title!("Analytical g2 with Ω/β = " * string(Ω/β))
    savefig("Figures/analytical_g2_" * string(Ω/β) * ".png")

end

### RUN CODE ###
β = 0.5
Ω = 2.3 * 0.5

end_time = 9 
time_steps = 1000

time_list, analytical_solution = analytical_array(end_time, time_steps, β, Ω)
plot_analytical_solution(time_list, analytical_solution, β, Ω)