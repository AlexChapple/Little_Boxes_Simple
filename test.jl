using Base.Threads

print("Number of threads: " * string(Threads.nthreads()),"\n")

# acc = Ref(0)

# @time @threads for i in 1:1000
#     acc[] += 1
# end

# print(acc)

# acc = Atomic{Int64}(0)

# @time @threads for i in 1:1000
#     atomic_add!(acc, 1)
# end

# print(acc)

###

function test()

    a = zeros(Threads.nthreads(),5) 

    @threads for i in 1:10000000
        a[Threads.threadid():Threads.threadid(), :] += [1 1 1 1 1]
    end
    
    b = zeros(1,5)
    for i in 1:Threads.nthreads()
        b += a[i:i, :]
    end

    print(b)

end

function test2()

    a = zeros(1,5)

    for i in 1:10000000
        a += [1 1 1 1 1]
    end

    print(a)


end

@time test()