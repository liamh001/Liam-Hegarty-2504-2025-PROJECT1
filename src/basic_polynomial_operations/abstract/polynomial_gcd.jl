#############################################################################
#############################################################################
#
# This file implements polynomial GCD for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for polynomials (of the same concrete subtype) modulo prime.
"""
function extended_euclid_alg(a::P, b::P, prime::Int) where {P <: Polynomial}
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(P), zero(P)
    old_t, t = zero(P), one(P)

    while !iszero(mod(r,prime))
        q = first(divide(old_r, r)(prime))
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t

    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end


"""
The extended euclid algorithm for polynomials (of the same concrete subtype) modulo prime.
"""
function square_free(f::P) where {P <: Polynomial}
    sq_fr_f = pseudo_quo(f, pseudo_gcd(f, derivative(f)))
    content = gcd(map(t -> t.coeff, sq_fr_f))
    return P( map(t -> Term(t.coeff ÷ content, t.degree), sq_fr_f) )
end

function prim_part(f::P) where {P <: Polynomial}
    iszero(f) && return f
    content = gcd(map(t -> t.coeff, f))
    return P( map(t -> Term(t.coeff ÷ content, t.degree), f) )
end

function pseudo_gcd(f::P, g::P) where {P <: Polynomial}
    while !(iszero(f) || iszero(g))
        f, g = g, pseudo_rem(f, g)
    end

    iszero(g) && return prim_part(f)
    return prim_part(g) # f = 0
end

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

"""
The GCD of two polynomials (of the same concrete subtype) modulo prime.
"""
gcd(a::P, b::P, prime::Int) where {P <: Polynomial} = extended_euclid_alg(a,b,prime) |> first