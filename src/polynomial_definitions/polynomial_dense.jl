#############################################################################
#############################################################################
#
# This file defines the `PolynomialDense` type with several operations 
#                                                                               
#############################################################################
#############################################################################

#########################################
# PolynomialDense type and construction #
#########################################

"""
A Polynomial type - designed to be for polynomials with integer coefficients.

This type utilises a dense representation for a polynomial. This means that
zero terms with degree less than the degree of the polynomial are stored.

E.g, for x^3 + 2x we store in memory:
    [Term(0, 0), Term(2, 1), Term(0, 2), Term(1, 3)]
"""
struct PolynomialDense <: Polynomial

    #A zero packed vector of terms
    #Terms are assumed to be in order with first term having degree 0, second degree 1, and so fourth
    #until the degree of the polynomial. The leading term (i.e. last) is assumed to be non-zero except 
    #for the zero polynomial where the vector is of length 1.
    #Note: at positions where the coefficient is 0, the power of the term is also 0 (this is how the Term type is designed)
    terms::Vector{Term}   
    
    #Inner constructor of 0 polynomial
    PolynomialDense() = new([zero(Term)])

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialDense(vt::Vector{Term})

        #Filter the vector so that there is not more than a single zero term
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term)]
        end

        max_degree = maximum((t)->t.degree, vt)
        terms = [zero(Term) for i in 0:max_degree] #First set all terms with zeros

        #now update based on the input terms
        for t in vt
            terms[t.degree + 1] = t #+1 accounts for 1-indexing
        end
        return new(terms)
    end
end

###########
# Display #
###########

# We don't need to override the `show` function from this section.
# The implementation from `Polynomial` is sufficient for our needs.`

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.
"""
iterate(p::PolynomialDense, state=1) = iterate(p.terms, state)

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialDense) = length(p.terms) 

"""
The leading term of the polynomial.
"""
leading(p::PolynomialDense)::Term = isempty(p.terms) ? zero(Term) : last(p.terms) 

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
# TODO/FIXME - THE BELOW IS LITERALLY A BUG
#Note that ideally this would throw and error if pushing another term of degree that is already in the polynomial
function push!(p::PolynomialDense, t::Term)
    if t.degree <= degree(p)
        p.terms[t.degree + 1] = t
    else
        append!(p.terms, zeros(Term, t.degree - degree(p)-1))
        push!(p.terms, t)
    end
    return p        
end

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialDense)::Term 
    popped_term = pop!(p.terms) #last element popped is leading coefficient

    while !isempty(p.terms) && iszero(last(p.terms))
        pop!(p.terms)
    end

    if isempty(p.terms)
        push!(p.terms, zero(Term))
    end

    return popped_term
end

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

# We again can rely on the Julia type system - the implementations for `Polynomial` will work for `PolynomialDense`

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same.

Note - even though this is done for `Polynomial`, we can override it for `PolynomialDense`
to leverage Julia's speed with vectors.
"""
==(p1::Polynomial, p2::Polynomial)::Bool = p1.terms == p2.terms

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

# Again, technically the abstract implementations for `Polynomial` will work correctly for `PolynomialDense`.
# However, this doesn't mean they are particularly efficient - if you wish you can override any particular
# operation and re-implement it for `PolynomialDense`.`