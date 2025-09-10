#############################################################################
#############################################################################
#
# This file implements operations for polynomials over a field, e.g. Zp[x].
#                                                                               
#############################################################################
#############################################################################

#= TASK LIST:
NOTE - the following functions will NOT work when C is an integer until you 
finish implementing all tasks (and even then, some functions such as gcd simply
cannot be implemented over the ring of integers - see bonus task 3 for details).

TODO (Task 2) 
    In this task, simply edit this file such that the first argument (::Type{C})
    to each function is also the type of the coefficients of the polynomials.
    I.e., we want the type signatures to contain:
        {C, D, P <: Polynomial{C, D}}

TODO (Task 5)
    In this task, override the unimplemented functions (after duplicating this 
    file as per the instructions) for C <: ZModP. When overriding, use the
    corresponding _mod_p function as a guide.

TODO (Task 6)
    Here, (in the duplicated file as per the instructions) override the factor
    function for C <: Integer and implement the Chinese Remainder Theorem (CRT).

=#

# TODO - ADD DOCSTRINGS

"""
Returns the factors of `f` in an array of tuples. 

For a given tuple (g, n) in the returned array, g is an irreducible factor
of f and the multiplicity of g in f is n.

The returned array may also contain a constant for reconstruction of the
leading coefficient.

NOTE: Override this in task 5 and 6.
"""
function factor(::Type{C}, f::P)::Vector{Tuple{P, Int}} where {C, P <: Polynomial}
    not_implemented_error(f, "factor")
end

"""
Returns the quotient and remainder of num divided by den (where num/den have the 
same concrete type). 

NOTE: Override this in task 5.
"""
function div_rem(::Type{C}, num::P, den::P)::Tuple{P, P} where {C, P <: Polynomial}
    not_implemented_error(num, "div_rem")
end

"""
Distinct degree factorization.

Given a square free polynomial `f` returns a list, `g` such that `g[k]` is a product 
of irreducible polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2, such that the 
product of the list is equal to `f`.
"""
function dd_factor(::Type{C}, f::P)::Array{P} where {C, P <: Polynomial}
    not_implemented_error(f, "dd_factor")
end

"""
Distinct degree split.

Returns a list of irreducible polynomials of degree `d` so that the product of
that list is the polynomial `f`.
"""
function dd_split(::Type{C}, f::P, d::Int)::Vector{P} where {C, P <: Polynomial}
    not_implemented_error(f, "dd_split")
end


""" 
Returns the quotient of num divided by den (where num/den have the 
same concrete type) 
"""
div(::Type{C}, num::P, den::P) where {C, P <: Polynomial} = first(div_rem(num, den))

""" 
Returns the remainder of num divided by den (where num/den have the 
same concrete type) 
"""
rem(::Type{C}, num::P, den::P) where {C, P <: Polynomial} = last(div_rem(num, den))

"""
The extended euclid algorithm for polynomials (of the same concrete subtype).
"""
function extended_euclid_alg(::Type{C}, f::P, g::P) where {C, P <: Polynomial}
    return ext_euclid_alg(f, g, rem, div)
end

"""
The greatest common divisor of two polynomials (of the same concrete subtype).
"""
gcd(::Type{C}, f::P, g::P) where {C, P <: Polynomial} = extended_euclid_alg(f, g) |> first

"""
Perfect field, yun's algorithm, square free
"""
function square_free(::Type{C}, f::P) where {C, P <: Polynomial}
    # TODO - See if I can implement this over any field using Yun's algorithm
    # TODO - Yun's algorithm holds over any field apparently
    # FIXME - TUTOR TO PROVIDE IMPLEMENTATION - MITTUN
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

"""
Compute the number of times g divides f.
"""
function multiplicity(::Type{C}, f::P, g::P)::Int where {C, P <: Polynomial}
    degree(gcd(f, g)) == 0 && return 0
    return 1 + multiplicity(div(f, g), g)
end
