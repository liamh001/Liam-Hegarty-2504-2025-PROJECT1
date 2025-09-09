#############################################################################
#############################################################################
#
# This file implements polynomial division for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""  Modular algorithm (f, g have the same concrete subtype).
f divide by g

f = q*g + r (mod p)

p is a prime

This must be overridden by concrete subtypes.
"""
function divide(num::P, den::P) where {P <: Polynomial}
    not_implemented_error(num, "divide")
end

"""
The quotient from polynomial division. Returns a function of an integer.
"""
÷(num::P, den::P) where {P <: Polynomial}  = (p::Int) -> first(divide(num,den)(p))

"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::P, den::P) where {P <: Polynomial} = (p::Int) -> last(divide(num,den)(p))

function pseudo_quo(f::P, g::P) where {P <: Polynomial}
    iszero(g) && error("Cannot divide by 0")
    m, n = degree(f), degree(g)
    lc_g_pow = one(P)
    ret_val = zero(P)
    lc_g = leading(g).coeff

    while !(m < n || iszero(f))
        x_pow = x_poly(P)^(m-n)
        lc_f = leading(f).coeff

        ret_val += lc_g_pow * lc_f * x_pow
        lc_g_pow *= lc_g
        f = (lc_g * f) - (lc_f * x_pow * g)

        m = degree(f)
    end

    ret_val += lc_g_pow * prim_part(f)
    return ret_val
end

function pseudo_rem(f::P, g::P) where {P <: Polynomial}
    m, n = degree(f), degree(g)
    while !(m < n || iszero(f))
        x_pow = x_poly(P)^(m-n)
        lc_f = leading(f).coeff
        lc_g = leading(g).coeff
        f = prim_part((lc_g * f) - (lc_f * x_pow * g))
        m = degree(f)
    end
    
    return f
end