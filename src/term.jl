#############################################################################
#############################################################################
#
# This file defines the Term type with several operations 
#                                                                               
#############################################################################
#############################################################################

##################################
# Term type and utility function #
##################################

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
        new(coeff, degree)
    end
end

######################
# Outer Constructors #
######################

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

"""
Creates the unit term.
"""
function one(::Type{Term{C, D}})::Term{C, D} where {C, D}
    Term(one(C), zero(D))
end
one(t::Term) = one(typeof(t))

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
function mod(t::Term{C, D}, p::Integer)::Term{C, D} where {C <: Integer, D}
    Term(mod(t.coeff, p), t.degree)
end

"""
Compute the derivative of a term.
"""
function derivative(t::Term{C, D})::Term{C, D} where {C, D}  
    Term{C, D}(t.coeff*C(t.degree),max(t.degree-one(D),zero(D)))
end

"""
Exact division when the coefficient type `C` is a field (e.g. Zp).

You will need to override this in Task 5 for ZModP specifically.
Do NOT modify this particular version of the function.
"""
function div(t1::Term{C, D}, t2::Term{C, D}) where {C, D} # Base.:÷
    not_yet_implemented_error(t1, "div")
end

"""
Divide two terms (with the same coefficient/degree types). 
Returns a function of an integer.

Note: You will need to override this for division where the coefficients are of type ZModP (Task 5)
There we can do exact division, so we can simply do `t1.coeff ÷ t2.coeff`.
"""
function div_mod_p(t1::Term{C, D}, t2::Term{C, D}, prime::Integer) where {C, D}
    @assert t1.degree ≥ t2.degree
    new_coeff = mod((mod(t1.coeff, prime) * int_inverse_mod(t2.coeff, prime)), prime)
    return Term(new_coeff, t1.degree - t2.degree)
end

"""
Integer divide a term by an integer.
"""
function div_mod_p(t::Term{C, D}, n::Integer, prime::Integer) where {C <: Integer, D} 
    return div_mod_p(t, Term(C(n), zero(D)), prime)
end

#############################
# Vectorization with a term #
#############################

""" Enable broadcasting an operation on a term """
broadcastable(t::Term) = Ref(t)




"""
Exact division for terms with ZModP coefficients.
"""
function div(t1::Term{ZModP{T, N}, D}, t2::Term{ZModP{T, N}, D}) where {T, N, D}
    @assert t1.degree ≥ t2.degree
    Term(t1.coeff ÷ t2.coeff, t1.degree - t2.degree)
end

"""
Refactored div_mod_p using ZModP conversion.
"""
function div_mod_p(t1::Term{C, D}, t2::Term{C, D}, prime::Integer) where {C <: Integer, D}
    @assert t1.degree ≥ t2.degree
    
    # Convert to ZModP terms
    t1_mod = Term(ZModP{C, prime}(t1.coeff), t1.degree)
    t2_mod = Term(ZModP{C, prime}(t2.coeff), t2.degree)
    
    # Perform division in Zp
    result_mod = div(t1_mod, t2_mod)
    
    # Convert back
    return Term(C(result_mod.coeff), result_mod.degree)
end