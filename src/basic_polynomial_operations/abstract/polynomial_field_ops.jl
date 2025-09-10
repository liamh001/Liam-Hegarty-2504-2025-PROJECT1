#############################################################################
#############################################################################
#
# This file implements operations for polynomials over a field, e.g. Zp[x].
#                                                                               
#############################################################################
#############################################################################

#= TASK LIST:

TODO (Task 2)

TODO (Task 5)

TODO (Task 6)

=#

# TODO - ADD FUNCTIONS ADD ABSTRACT LEVEL FOR F[X]
# TODO - ADD DOCSTRINGS

function div_rem(num::P, den::P)::Tuple{P, P} where {P <: Polynomial}
    not_implemented_error(num, "div_rem")
end

div(num::P, den::P) where {P <: Polynomial} = first(div_rem(num, den))
rem(num::P, den::P) where {P <: Polynomial} = last(div_rem(num, den))

function extended_euclid_alg(f::P, g::P) where {P <: Polynomial}
    return ext_euclid_alg(f, g, rem, div)
end

gcd(a::P, b::P) where {P <: Polynomial} = extended_euclid_alg(a, b) |> first

# TODO - IMPLEMENT PROPERLY OR LEAVE UNIMPLEMENTED!
function square_free(f::P) where {P <: Polynomial}
    not_implemented_error(f, "square_free")
    # fmod_p = mod(f, prime)

    # min_deg = last(fmod_p).degree
    # vt = filter(t -> !iszero(t), collect(f))
    # fmod_p = P( map(t -> Term(t.coeff, t.degree - min_deg), vt) )

    # # Now compute the gcd modulo a prime
    # der_fmod_p = mod(derivative(fmod_p), prime)
    # gcd_f_der_f = gcd_mod_p(fmod_p, der_fmod_p, prime)

    # iszero(gcd_f_der_f) && return fmod_p * (min_deg > zero(min_deg) ? x_poly(P) : one(P))

    # sqr_free = div_mod_p(fmod_p, gcd_f_der_f, prime)

    # # Add the factor of `x` back in if there was one
    # if min_deg > zero(min_deg) 
    #     sqr_free *= x_poly(P)
    # end

    # return sqr_free
end