#############################################################################
#############################################################################
#
# This file defines the `PolynomialSparse` type with several operations 
#                                                                               
#############################################################################
#############################################################################

#########################################
# PolynomialSparse type and construction #
#########################################

"""
A Polynomial type - designed to be for polynomials with sparse representation.

This type utilises a heap-based sparse representation for a polynomial. This means that
only non-zero terms are stored, ordered by degree in a max-heap.

E.g, for x^100 + 2x^50 + 3 we only store:
    Heap([Term(1, 100), Term(2, 50), Term(3, 0)])
"""
struct PolynomialSparse{C, D} <: Polynomial{C, D}
    terms::Heap{Term{C, D}}
    
    function PolynomialSparse{C,D}() where {C,D}
        new(Heap{Term{C, D}}())
    end

    function PolynomialSparse{C, D}(vt::Vector{Term{C, D}}) where {C, D}
        non_zero_terms = filter(t -> !iszero(t), vt)
        
        combined_terms = Dict{D, C}()
        for t in non_zero_terms
            if haskey(combined_terms, t.degree)
                combined_terms[t.degree] += t.coeff
            else
                combined_terms[t.degree] = t.coeff
            end
        end
        
        final_terms = Term{C, D}[]
        for (deg, coeff) in combined_terms
            if !iszero(coeff)
                push!(final_terms, Term(coeff, deg))
            end
        end
        
        h = Heap(final_terms)
        new(h)
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

function iterate(p::PolynomialSparse{C,D}, state=nothing) where {C,D}
    if state === nothing
        sorted_terms = sort([t for t in p.terms.data if !iszero(t)], by = t -> t.degree)
        if isempty(sorted_terms)
            return nothing
        end
        return (sorted_terms[1], (sorted_terms, 2))
    else
        sorted_terms, idx = state
        if idx > length(sorted_terms)
            return nothing
        end
        return (sorted_terms[idx], (sorted_terms, idx + 1))
    end
end

##############################
# Queries about a polynomial #
##############################

function length(p::PolynomialSparse)
    count(!iszero, p.terms.data)
end

function leading(p::PolynomialSparse{C,D})::Term{C,D} where {C,D}
    isempty(p.terms) ? zero(Term{C, D}) : peek(p.terms)
end

function last(p::PolynomialSparse{C,D}) where {C,D}
    iszero(p) && return zero(Term{C, D})
    min_term = nothing
    for t in p.terms.data
        if !iszero(t) && (min_term === nothing || t.degree < min_term.degree)
            min_term = t
        end
    end
    return min_term === nothing ? zero(Term{C, D}) : min_term
end

################################
# Pushing and popping of terms #
################################

function push!(p::PolynomialSparse{C,D}, t::Term{C,D}) where {C,D}
    if iszero(t)
        return p
    end
    
    for i in 1:length(p.terms.data)
        if p.terms.data[i].degree == t.degree
            error("Cannot push a term with degree $(t.degree) that already exists")
        end
    end
    
    push!(p.terms, t)
    return p
end

function pop!(p::PolynomialSparse{C,D})::Term{C,D} where {C,D}
    isempty(p.terms) ? zero(Term{C, D}) : pop!(p.terms)
end

#################################
# Queries about two polynomials #
#################################

function ==(p1::PolynomialSparse{C,D}, p2::PolynomialSparse{C,D})::Bool where {C,D}
    terms1 = Dict{D, C}()
    terms2 = Dict{D, C}()
    
    for t in p1.terms.data
        if !iszero(t)
            terms1[t.degree] = t.coeff
        end
    end
    
    for t in p2.terms.data
        if !iszero(t)
            terms2[t.degree] = t.coeff
        end
    end
    
    return terms1 == terms2
end

#################################################################
# Special operations for sparse polynomial efficiency           #
#################################################################

function iszero(p::PolynomialSparse)::Bool
    for t in p.terms.data
        if !iszero(t)
            return false
        end
    end
    return true
end

function degree(p::PolynomialSparse{C,D})::D where {C,D}
    iszero(p) ? zero(D) : leading(p).degree
end