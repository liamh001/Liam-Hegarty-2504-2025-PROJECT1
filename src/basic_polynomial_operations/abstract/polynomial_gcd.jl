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
    if iszero(g) 
        return prim_part(f)
    elseif iszero(f)
        return prim_part(g)
    end

    r = pseudo_rem(f, g)

    return pseudo_gcd(g, r)
end

function pseudo_quo(f::P, g::P) where {P <: Polynomial}
    iszero(g) && error("Cannot divide by 0")
    m, n = degree(f), degree(g)
    (m < n || iszero(f)) && return prim_part(f)

    lc_f = leading(f).coeff
    lc_g = leading(g).coeff
    t = m - n
    x_pow = x_poly(P)^t

    f1 = (lc_g * f) - (lc_f * x_pow * g)

    return (lc_g * pseudo_quo(f1, g)) + (lc_f * x_pow) 
end

function pseudo_rem(f::P, g::P) where {P <: Polynomial}
    m, n = degree(f), degree(g)
    (m < n || iszero(f)) && return prim_part(f)

    lc_f = leading(f).coeff
    lc_g = leading(g).coeff
    t = m - n
    x_pow = x_poly(P)^t

    f1 = (lc_g * f) - (lc_f * x_pow * g)

    return pseudo_rem(f1, g)
end

"""
The GCD of two polynomials (of the same concrete subtype) modulo prime.
"""
gcd(a::P, b::P, prime::Int) where {P <: Polynomial} = extended_euclid_alg(a,b,prime) |> first