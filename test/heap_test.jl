
function test_heap(;N::Int = 10^4)
    # Test empty constructor
    h = Heap{Int}()
    v = rand(1:N, N)
    for i in 1:N
        push!(h, v[i])
    end
    sort!(v, rev=true)

    @assert peek(h) == maximum(v)
    @assert peek(h) == maximum(v) # Check peek does not affect the heap
    heap_sorted = popall!(h)
    !isempty(h) && error("After popping all elements the heap should be empty")
    v != heap_sorted && error("Popping all from the heap should return a sorted vector")
    # println("test_heap - push!/popall!/empty constructor - PASSED")

    # push!/pop!/peek/isempty
    x = 1
    push!(h, x)
    x = pop!(h)
    push!(h, x)
    @assert x == 1
    @assert peek(h) == 1
    pop!(h)
    @assert isempty(h)
    # println("test_heap - push!/pop!/peek/isempty - PASSED")

    # non-destructive safe vector constructor
    v = rand(1:N, N)
    h = Heap(v)
    @assert peek(h) == maximum(v)
    heap_sorted = popall!(h)
    sort!(v, rev=true)
    v != heap_sorted && error("Popping all from heap should return vector in sorted order")
    @assert isempty(h)
    # println("test_heap - non-destructive safe vector constructor - PASSED")

    # destructive safe vector constructor
    v = rand(1:N, N)
    v_cpy = sort(v, rev=true)
    h = Heap!(v)
    @assert peek(h) == maximum(v)
    heap_sorted = popall!(h)
    @assert v_cpy == heap_sorted
    @assert isempty(h)
    @assert isempty(v)
    # println("test_heap - destructive safe vector constructor - PASSED")

    # unsafe vector constructor
    v = rand(1:N, N)
    h = unsafe_vec_to_heap(v)
    @assert peek(h) == v[1]
    @assert length(h) == N
    # println("test_heap - unsafe vector constructor - PASSED")

    println("test_heap - PASSED")
end
