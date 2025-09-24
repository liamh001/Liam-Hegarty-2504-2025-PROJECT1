#############################################################################
#############################################################################
#
# This file implements polynomial addition for sparse polynomials.
#                                                                               
#############################################################################
#############################################################################

function +(p::PolynomialSparse{C, D}, t::Term{C, D}) where {C, D}
    iszero(t) && return deepcopy(p)
    
    p_new = PolynomialSparse{C, D}()
    found = false
    
    for term in p.terms.data
        if !iszero(term)
            if term.degree == t.degree
                new_coeff = term.coeff + t.coeff
                if !iszero(new_coeff)
                    push!(p_new.terms, Term(new_coeff, t.degree))
                end
                found = true
            else
                push!(p_new.terms, term)
            end
        end
    end
    
    if !found
        push!(p_new.terms, t)
    end
    
    return p_new
end

function +(p1::PolynomialSparse{C,D}, p2::PolynomialSparse{C,D})::PolynomialSparse{C,D} where {C,D}
    terms_dict = Dict{D, C}()
    
    for t in p1.terms.data
        if !iszero(t)
            terms_dict[t.degree] = t.coeff
        end
    end
    
    for t in p2.terms.data
        if !iszero(t)
            if haskey(terms_dict, t.degree)
                terms_dict[t.degree] += t.coeff
            else
                terms_dict[t.degree] = t.coeff
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