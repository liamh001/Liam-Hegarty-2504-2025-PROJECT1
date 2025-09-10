#############################################################################
#############################################################################
#
# This file implements pseudo polynomial division for abstract polynomials
# This file is ONLY here for bonus task 3 - do not attempt to use these
# functions otherwise.
#                                                                               
#############################################################################
#############################################################################

# """
# A square free (and primitive) polynomial over Z[x] via pseudo quotient and remainder..
# This function utilises pseudo-division over the ring Z[x] and may be
# subject to coefficient swell leading to extremely bad performance.

# Details at:
# https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Pseudo-remainder_sequences 

# E.g.,
#     f = (x + a1)^e1 * (x + a2)^e2 * ... * (x + an)^en
#     square_free(f) = (x + a1) * (x + a2) * ... * (x + an)

# Do NOT use this for `Zp[x]`.
# """
# function pseudo_square_free(f::P) where {P <: Polynomial}
#     sq_fr_f = pseudo_quo(f, pseudo_gcd(f, derivative(f)))
#     return prim_part(sq_fr_f)
# end


# """
# Computes the pseudo gcd of two polynomials (of the same concrete subtype).

# We cannot compute the gcd over Z[x] as Z is not a field. We can either compute the gcd in the related fraction
# field (i.e., Q[x]) or use the `pseudo_gcd`. The gcd may have fractional coefficients, to go from this to the
# pseudo_gcd we multiply by the lcm of the denominators of these fractions. This however, can still be slow
# as the numerator and denominator can respectively become very large.

# The following algorithm also computes the pseudo_gcd without resorting to Q[x] (i.e., we stay in Z[x]).
# """
# function pseudo_gcd(f::P, g::P) where {P <: Polynomial}
#     while !(iszero(f) || iszero(g))
#         f, g = g, pseudo_rem(f, g)
#     end

#     iszero(g) && return prim_part(f)
#     return prim_part(g) # f = 0
# end
