using PyPlot
using Random
using LinearAlgebra
pygui(:qt)


a = LinRange(0,10,10)
b = LinRange(0,10,20)

print(size(a), size(b))

c = zeros(size(a)[1], size(b)[1])

for i in 1:size(a)[1]
    for j in 1:size(b)[1]
        c[i,j] = i^2
    end
end

b = reverse(b)
c = c[end:-1:1, :]


surf(a,b,c')
xlabel("Γ")
ylabel("γt") 
zlabel("w(γt)")
plt.gca().invert_yaxis()
# plt.savefig("Figures/testing_surf.png", dpi=600)
