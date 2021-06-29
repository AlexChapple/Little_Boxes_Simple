using LinearAlgebra

#= 
Demonstrates how binning works in the waiting distribution code. 
=#

a = 8

b = LinRange(0,9,10)
c = zeros(size(b))

if a < b[2]
    c[1] += 1
    print("first clause executed")
elseif a >= b[end-1] && a <= b[end]
    c[end] += 1
    print("second clause executed")
else
    for k in 2:(size(b)[1]-1)
        if a >= b[k] && a < b[k+1]
            c[k] += 1
            print("third clause executed")
            break
        end
    end
end

for i in b
    print(i, "  ")
end
print("\n")
for i in c
    print(i, "  ")
end





