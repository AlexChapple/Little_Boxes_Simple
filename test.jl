using Base.Threads

print(Threads.nthreads())

acc = Ref(0)

@time @threads for i in 1:1000
    acc[] += 1
end

print(acc)

acc = Atomic{Int64}(0)

@time @threads for i in 1:1000
    atomic_add!(acc, 1)
end

print(acc)