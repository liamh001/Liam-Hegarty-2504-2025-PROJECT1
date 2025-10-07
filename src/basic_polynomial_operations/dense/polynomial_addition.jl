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
function +(p::PolynomialDense{C, D}, t::Term{C, D}) where {C, D}
    p = deepcopy(p)
    if t.degree > degree(p)
        push!(p, t)
    else
        new_coeff = p.terms[t.degree + 1].coeff + t.coeff
        p.terms[t.degree + 1] = Term(new_coeff, t.degree)
    end
    return trim!(p)
end

"""
Add two dense polynomials efficiently by working directly with coefficient arrays.
"""
function +(p1::PolynomialDense{C,D}, p2::PolynomialDense{C,D})::PolynomialDense{C,D} where {C,D}
    deg1 = degree(p1)
    deg2 = degree(p2)
    max_deg = max(deg1, deg2)
    
    result_terms = [Term(zero(C), D(i)) for i in 0:max_deg]
    
    for i in 1:min(length(p1.terms), deg1+1)
        result_terms[i] = Term(result_terms[i].coeff + p1.terms[i].coeff, D(i-1))
    end
    
    for i in 1:min(length(p2.terms), deg2+1)
        result_terms[i] = Term(result_terms[i].coeff + p2.terms[i].coeff, D(i-1))
    end
    
    return trim!(PolynomialDense{C,D}(result_terms))
end


