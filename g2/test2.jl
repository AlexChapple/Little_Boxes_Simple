using LinearAlgebra

#= 
Demonstrates how binning works in the waiting distribution code. 
=#

d = [1 2 3 4 5 6 7 8]


b = LinRange(0,9,5)
c = zeros(size(b))

for a in d

    if a < b[2]
        c[1] += 1
        # print("first clause executed")
    elseif a >= b[end-1] && a <= b[end]
        c[end-1] += 1
        # print("second clause executed")
    else
        for k in 2:(size(b)[1]-1)
            if a >= b[k] && a < b[k+1]
                c[k] += 1
                # print("third clause executed")
                break
            end
        end
    end
end

for i in b
    print(i)
end

plot(b,c)






