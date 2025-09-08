#############################################################################
#############################################################################
#
# This file implements polynomial addition for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.
"""
function +(p::Polynomial, t::Term)
    not_implemented_error(p, "Polynomial + Term")
end
+(t::Term, p::Polynomial) = p + t

"""
Add two polynomials.
"""
function +(p1::Polynomial, p2::Polynomial)::Polynomial
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

"""
Add a polynomial and an integer.
"""
+(p::Polynomial, n::Int) = p + Term(n,0)
+(n::Int, p::Polynomial) = p + Term(n,0)