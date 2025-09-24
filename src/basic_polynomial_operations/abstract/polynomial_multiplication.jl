#############################################################################
#############################################################################
#
# This file implements polynomial multiplication for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two polynomials (of the same concrete subtype).
"""
function *(p1::P, p2::P)::P where {P <: Polynomial}
    p_out = P()
    for t in p1
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end

"""
Power of a polynomial.
"""
function ^(p::P, n::Integer)::P where {P <: Polynomial}
    n < 0 && error("No negative power")
    n == 0 && return one(p)
    n == 1 && return deepcopy(p)
    
    # Special cases for O(1) operations
    iszero(p) && return zero(p)  # 0^n = 0
    p == one(p) && return one(p)  # 1^n = 1
    
    # Check if p is just x (single term with coefficient 1 and degree 1)
    if length(p) == 1 && leading(p).coeff == 1 && degree(p) == 1
        # x^n can be constructed directly
        return P([Term(one(leading(p).coeff), n * degree(p))])
    end
    
    # Repeated squaring for general case
    result = one(p)
    base = deepcopy(p)
    
    while n > 0
        if n % 2 == 1
            result *= base
        end
        base *= base
        n ÷= 2
    end
    
    return result
end

