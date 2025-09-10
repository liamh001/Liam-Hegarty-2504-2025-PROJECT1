#############################################################################
#############################################################################
#
# This file implements polynomial division for dense polynomials.
#                                                                               
#############################################################################
#############################################################################

# TODO - MAKE THIS ABSTRACT
"""  Modular algorithm.
f divide by g

f = q*g + r

p is a prime
"""
function divide_mod_p(num::PolynomialDense, den::PolynomialDense, prime::Int)
    f, g = mod(num,prime), mod(den,prime)
    degree(f) < degree(num) && return nothing 
    iszero(g) && throw(DivideError())
    q = PolynomialDense()
    prev_degree = degree(f)
    while degree(f) ≥ degree(g) 
        h = PolynomialDense( (leading(f) ÷ leading(g))(prime) )  #syzergy 
        f = mod((f - h*g), prime)
        q = mod((q + h), prime)  
        prev_degree == degree(f) && break
        prev_degree = degree(f)
    end
    @assert iszero( mod((num  - (q*g + f)),prime))
    return q, f
end


# We won't re-implement any of these functions for dense polynomials, the abstract versions will 
# produce the correct result.

# ÷(num::Polynomial, den::Polynomial)
# rem_mod_p(num::Polynomial, den::Polynomial, prime::Int)