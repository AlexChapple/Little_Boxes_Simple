using PyPlot; pyplot();
using Random
pygui(:qt)


a = [1,2,3,4,5]
b = [1,2,3,4,5]

c = zeros(size(a)[1], size(b)[1])


for i in 1:size(a)[1]
    for j in 1:size(b)[1]
        c[i,j] = rand()
    end
end


surface(a,b,c)
