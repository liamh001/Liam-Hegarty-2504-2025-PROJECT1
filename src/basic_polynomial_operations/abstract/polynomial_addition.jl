#############################################################################
#############################################################################
#
# This file implements polynomial addition for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.

This must be overridden by concrete subtypes.
"""
function +(p::Polynomial, t::Term)
    not_implemented_error(p, "Polynomial + Term")
end
+(t::Term, p::Polynomial) = p + t

"""
Add two polynomials of the same concrete subtype.

Note: This operation may be slow for some concrete subtypes. You may wish to override this to factor
in the details of your polynomial representation when implementing your concrete subtype.
"""
function +(p1::P, p2::P)::P where {P <: Polynomial}
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

"""
Add a polynomial and an integer.
"""
+(p::Polynomial, n::Integer) = p + Term(n,0)
+(n::Integer, p::Polynomial) = p + Term(n,0)