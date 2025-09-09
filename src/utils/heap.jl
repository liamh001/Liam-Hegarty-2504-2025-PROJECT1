#############################################################################
#############################################################################
#
# This file implements the Heap data structure 
#                                                                               
#############################################################################
#############################################################################


""" 
A struct to represent a basic array-based max-heap.

T -> can be any type which implements comparison (e.g., integers, reals etc.)

The underlying array representation can be accessed via h.data where h::Heap.
"""
struct Heap{T}
    data::Vector{T}

    """ 
    Inner constructor to create a heap from an array.

    WARNING:
        Does NOT check whether the heap property is maintained and does not enforce this.
        You must ensure the heap property is satisfied (e.g., the array is sorted in descending order)
        before using this constructor.  
    """
    Heap{T}(v::Vector{T}) where {T} = new(v) 
end

#####################
# Heap constructors #
#####################

"""
Given a vector, turns this into a heap (maintaining heap order)
"""
function vec_to_heap(v::Vector{T}) where {T}
    h = Heap{T}(v)
    _heapify!(h)
    return h
end

""" Create an empty heap """
Heap{T}() where {T} = Heap(T[])

""" Convenience function to create a heap and infer the type """
Heap(v::Vector{T}) where {T} = Heap{T}(v)

###########
# Display #
###########

show(io, h::Heap) = print(io, "Heap(", h.data, ")")

#####################
# Querying the heap #
#####################

""" Get the number of elements in the heap. """
length(h::Heap) = length(h.data)

""" 
Push a new element onto the heap.
The heap will automatically handle maintaining heap order
"""
function push!(h::Heap{T}, x::T) where {T}
    push!(h.data, x)
    _up_heap(h, length(h.data))
end

""" 
Pop (remove and return) the maximum element from the heap.
The heap will automatically handle maintaining heap order.

Note: This function will error if called when the heap is empty - check first with isempty(h)
"""
function pop!(h::Heap)
    isempty(h) && error("Heap is empty, cannot pop!")
    
    # Swap first/last elements
    max_val = h.data[1]
    h.data[1] = pop!(h.data)
    
    # Restore heap order
    if !isempty(h)
        _down_heap!(h, 1)
    end
    
    return max_val
end

"""
Returns the maximum element from the heap without removing it.

Note: This function will error if called when the heap is empty - check first with isempty(h)
"""
function peek(h::Heap) 
    isempty(h) && error("Heap is empty, cannot peek")
    return h.data[1]
end

####################################################
# Internal helper functions to maintain heap order #
####################################################

"""
Reorders a vector in place to satisfy the heap property.
"""
function _heapify!(h::Heap)
    n = length(h.data)
    for i in div(n, 2):-1:1
        _down_heap!(h, i)
    end
end

""" Internal helper function to restore heap order upwards """
function _up_heap!(h::Heap, i::Integer)
    parent_index = div(i, 2)
    while i > 1 && h.data[i] < h.data[parent_index]
        h.data[i], h.data[parent_index] = h.data[parent_index], h.data[i]
        i = parent_index
        parent_index = div(i, 2)
    end
end

""" Internal helper function to restore heap order downwards """
function _down_heap!(h::Heap, i::Integer)
    n = length(h.data)
    while (left_child = 2 * i) <= n
        j = left_child
        if (right_child = 2 * i + 1) <= n && h.data[right_child] < h.data[left_child]
            j = right_child
        end

        h.data[i] < h.data[j] && break
        h.data[i], h.data[j] = h.data[j], h.data[i]
        i = j
    end
end
