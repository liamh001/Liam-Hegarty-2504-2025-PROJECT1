#############################################################################
#############################################################################
#
# This file defines the Term type with several operations 
#                                                                               
#############################################################################
#############################################################################

##############################
# Term type and construction #
##############################

"""
A term.

For Term{C, D},
    C: The type of the coefficient of the term.
    D: The type of the degree of the term.
"""
struct Term{C <: Integer, D <: Integer} #structs are immutable by default
    coeff::C
    degree::D

    function Term{C, D}(coeff::C, degree::D) where {C, D}
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new(coeff, degree) : new(zero(C), zero(D))
    end
end

""" Convenience outer constructor to infer coefficient/degree types. """
function Term(coeff::C, degree::D) where {C, D}
    Term{C, D}(coeff, degree)
end

"""
Creates the zero term.
"""
function zero(::Type{Term{C, D}})::Term{C, D} where {C, D} 
    Term(zero(C), zero(D))
end
zero(t::Term) = zero(typeof(t))
zero(::Type{Term}) = zero(Term{Int, Int}) # Default constructor

"""
Creates the unit term.
"""
function one(::Type{Term{C, D}})::Term{C, D} where {C, D}
    Term(one(C), zero(D))
end
one(t::Term) = one(typeof(t))
one(::Type{Term}) = one(Term{Int, Int}) # Default constructor

###########
# Display #
###########

"""
Show a term.
"""
show(io::IO, t::Term) = print(io, "$(t.coeff)⋅x^$(t.degree)")

########################
# Queries about a term #
########################

"""
Check if a term is 0.
"""
iszero(t::Term)::Bool = iszero(t.coeff)

"""
Compare two terms (with the same coefficient/degree types).
"""
function isless(t1::Term{C, D}, t2::Term{C, D})::Bool  where {C, D}
    t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  
end

function (==)(t1::Term{C, D}, t2::Term{C, D}) where {C, D}
    t1.coeff == t2.coeff && t1.degree == t2.degree
end

"""
Evaluate a term at a point x.
"""
function evaluate(t::Term, x::S) where {S <: Number} 
    t.coeff * x^t.degree
end

##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree (with the same coefficient/degree types).
"""
function +(t1::Term{C, D},t2::Term{C, D})::Term{C, D} where {C, D}
    @assert t1.degree == t2.degree
    Term(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
function -(t::Term{C, D},)::Term{C, D} where {C, D} 
    Term(-t.coeff,t.degree)  
end

"""
Subtract two terms with the same degree (and the same coefficient/degree types).
"""
function -(t1::Term{C, D}, t2::Term{C, D})::Term{C, D} where {C, D}
    t1 + (-t2) 
end

"""
Multiply two terms (with the same coefficient/degree types).
"""
function *(t1::Term{C, D}, t2::Term{C, D})::Term{C, D} where {C, D}
    Term(t1.coeff * t2.coeff, t1.degree + t2.degree)
end

"""
Multiply a term by a constant.
"""
function *(t::Term{C, D}, n::S)::Term{C, D} where {C, D, S <: Integer}
    Term(t.coeff * n, t.degree)
end
function *(n::S, t::Term{C, D})::Term{C, D} where {C, D, S <: Integer}
    t * n
end

"""
Power of a term.
"""
function ^(t::Term{C, D}, n::S)::Term{C, D} where {C, D, S <: Integer}
    Term(t.coeff^n, t.degree*D(n))
end

"""
Compute the symmetric mod of a term with an integer.
"""
function mod(t::Term{C, D}, p::Int)::Term{C, D} where {C <: Integer, D}
    Term(mod(t.coeff,p), t.degree)
end

"""
Compute the derivative of a term.
"""
function derivative(t::Term{C, D})::Term{C, D} where {C, D}  
    Term{C, D}(t.coeff*C(t.degree),max(t.degree-one(D),zero(D)))
end

"""
Divide two terms (with the same coefficient/degree types). 
Returns a function of an integer.

Note: You will need to override this for division where the coefficients are of type ZModP (Task 5)
There we can do exact division, so we can simply do `t1.coeff ÷ t2.coeff`.
"""
function ÷(t1::Term{C, D}, t2::Term{C, D}) where {C, D} #\div + [TAB]
    @assert t1.degree ≥ t2.degree
    f(p::Int)::Term{C, D} = Term{C, D}(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p), t1.degree - t2.degree)
end

"""
Integer divide a term by an integer.
"""
function ÷(t::Term{C, D}, n::C) where {C <: Integer, D} 
    t ÷ Term(C(n), zero(D))
end

#############################
# Vectorisation with a term #
#############################

""" Enable broadcasting an operation on a term """
broadcastable(t::Term) = Ref(t)
