#############################################################################
#############################################################################
#
# This file contains defines unit tests for heap operations
#                                                                               
#############################################################################
#############################################################################


"""
Tests heap interface.
"""
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

    # non-destructive pop all
    v = rand(1:N, N)
    h = Heap(v)
    heap_sorted = popall(h)
    sort!(v, rev=true)
    v != heap_sorted && error("Popping all from heap should return vector in sorted order")
    @assert length(h) == N
    @assert first(v) == peek(h)
    @assert !isempty(h)
    println("test_heap - popall non-destructive vector constructor - PASSED")

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

    # Non-destructive heap mapping
    v = rand(1:N, N)
    h1 = Heap(v)
    h2 = map_heap(h1, t -> 2 * t)
    @assert 2 * peek(h1) == peek(h2)
    v1 = popall!(h1)
    v2 = popall!(h2)
    @assert [2*t for t in v1] == v2
    @assert isempty(h1) && isempty(h2)
    # println("test_heap - non-destructive map heap - PASSED")

    # Destructive heap mapping
    v = rand(1:N, N)
    h = Heap(v)
    h_cpy = deepcopy(h)
    map_heap!(h, t -> 2 * t)
    @assert length(h) == length(h_cpy) == N
    @assert 2 * peek(h_cpy) == peek(h)
    v1 = popall!(h)
    v2 = popall!(h_cpy)
    @assert [2*t for t in v2] == v1
    @assert isempty(h1) && isempty(h2)
    # println("test_heap - destructive map heap - PASSED")

    println("test_heap - PASSED")
end
