#############################################################################
#############################################################################
#
# This file implements polynomial addition for dense polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.
"""
function +(p::PolynomialDense, t::Term)
    p = deepcopy(p)
    if t.degree > degree(p)
        push!(p, t)
    else
        if !iszero(p.terms[t.degree + 1]) #+1 is due to indexing
            p.terms[t.degree + 1] += t
        else
            p.terms[t.degree + 1] = t
        end
    end

    return trim!(p)
end

# FIXME
# """
# Add two polynomials. 
# """
# function +(p1::PolynomialDense, p2::PolynomialDense)::PolynomialDense
#     min_degree = min(degree(p1), degree(p2))
#     lower_order_terms = p1.terms[1:min_degree] .+ p2.terms[1:min_degree]

#     if degree(p1) == degree(p2)
#         return PolynomialDense(lower_order_terms)
#     end

#     if degree(p1) > degree(p2)
#         higher_order_terms = p1.terms[min_degree+1:end]
#     else degree(p1) < degree(p2)
#         higher_order_terms = p2.terms[min_degree+1:end]
#     end

#     return PolynomialDense([higher_order_terms; lower_order_terms])
# end
