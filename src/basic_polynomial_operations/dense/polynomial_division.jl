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

Returns quotient and remainder (q, r)
"""
function div_rem_mod_p(num::PolynomialDense, den::PolynomialDense, prime::Int)::Tuple{PolynomialDense, PolynomialDense}
    f, g = mod(num,prime), mod(den,prime)
    degree(f) < degree(num) && return nothing # FIXME/TODO - This breaks factor I think - why would we return nothing in place of a Tuple? - Mittun
    iszero(g) && throw(DivideError())
    q = PolynomialDense()
    prev_degree = degree(f)
    while degree(f) ≥ degree(g) 
        h = PolynomialDense( div_mod_p(leading(f), leading(g), prime) )  #syzergy 
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

# div_mod_p(num::Polynomial, den::Polynomial, prime::Int)
# rem_mod_p(num::Polynomial, den::Polynomial, prime::Int)