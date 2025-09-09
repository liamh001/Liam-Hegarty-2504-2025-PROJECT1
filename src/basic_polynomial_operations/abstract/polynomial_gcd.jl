#############################################################################
#############################################################################
#
# This file implements polynomial GCD for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Computes the pseudo gcd of two polynomials (of the same concrete subtype).

We cannot compute the gcd over Z[x] as Z is not a field. We can either compute the gcd in the related fraction
field (i.e., Q[x]) or use the `pseudo_gcd`. The gcd may have fractional coefficients, to go from this to the
pseudo_gcd we multiply by the lcm of the denominators of these fractions.

The following algorithm also computes the pseudo_gcd without resorting to Q[x] (i.e., we stay in Z[x]).
"""
function pseudo_gcd(f::P, g::P) where {P <: Polynomial}
    while !(iszero(f) || iszero(g))
        f, g = g, pseudo_rem(f, g)
    end

    iszero(g) && return prim_part(f)
    return prim_part(g) # f = 0
end

"""
The extended euclid algorithm for polynomials (of the same concrete subtype) modulo prime.

Note, when working with polynomials over Z[x] we cannot compute the extended euclidean algorithm (EEA) directly.
This is because we do not have division in Z[x]. Hence, we must drop to the field Zp[x] and compute the EEA
with respect to some prime. 

When you implement polynomials over Zp you can create your own version of this function that does not require 
a prime as input.
"""
function extended_euclid_alg_mod_p(a::P, b::P, prime::Int) where {P <: Polynomial}
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
The GCD of two polynomials (of the same concrete subtype) modulo prime.

Again, when you implement polynomials over Zp you can create your own version of this function that does not
require a prime as input.
"""
gcd_mod_p(a::P, b::P, prime::Int) where {P <: Polynomial} = extended_euclid_alg_mod_p(a,b,prime) |> first
