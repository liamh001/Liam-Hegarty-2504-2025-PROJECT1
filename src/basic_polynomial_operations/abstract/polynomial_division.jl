#############################################################################
#############################################################################
#
# This file implements polynomial division for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""  Modular algorithm.
f divide by g

f = q*g + r

p is a prime

Returns quotient and remainder (q, r)
"""
function div_rem_mod_p(num::P, den::P, prime::Int)::Tuple{P, P} where {P <: Polynomial}
    f, g = mod(num,prime), mod(den,prime)
    degree(f) < degree(num) && return nothing # FIXME - Mittun
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

"""
The quotient from polynomial division modulo a prime. 
"""
div_mod_p(num::P, den::P, prime::Int) where {P <: Polynomial} = first(div_rem_mod_p(num, den, prime))

"""
The remainder from polynomial division modulo a prime.
"""
rem_mod_p(num::P, den::P, prime::Int) where {P <: Polynomial} = last(div_rem_mod_p(num, den, prime))
