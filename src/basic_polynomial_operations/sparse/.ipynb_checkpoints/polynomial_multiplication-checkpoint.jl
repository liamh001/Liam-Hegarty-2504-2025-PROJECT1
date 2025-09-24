

#############################################################################
#############################################################################
#
# This file implements polynomial multiplication for sparse polynomials.
#                                                                               
#############################################################################
#############################################################################

function *(t::Term{C,D}, p::PolynomialSparse{C,D})::PolynomialSparse{C,D} where {C,D}
    iszero(t) && return PolynomialSparse{C,D}()
    
    result_terms = Term{C,D}[]
    for term in p.terms.data
        if !iszero(term)
            push!(result_terms, t * term)
        end
    end
    
    return PolynomialSparse{C,D}(result_terms)
end
*(p::PolynomialSparse, t::Term) = t * p

function *(p1::PolynomialSparse{C,D}, p2::PolynomialSparse{C,D})::PolynomialSparse{C,D} where {C,D}
    terms_dict = Dict{D, C}()
    
    for t1 in p1.terms.data
        if !iszero(t1)
            for t2 in p2.terms.data
                if !iszero(t2)
                    product_term = t1 * t2
                    deg = product_term.degree
                    if haskey(terms_dict, deg)
                        terms_dict[deg] += product_term.coeff
                    else
                        terms_dict[deg] = product_term.coeff
                    end
                end
            end
        end
    end
    
    result_terms = Term{C, D}[]
    for (deg, coeff) in terms_dict
        if !iszero(coeff)
            push!(result_terms, Term(coeff, deg))
        end
    end
    
    return PolynomialSparse{C, D}(result_terms)
end