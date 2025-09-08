#############################################################################
#############################################################################
#
# This file defines the abstract polynomial type alongside a number of common 
# operations that concrete subtypes should implement.
# Some functionality can be implemented at this abstract level.
# However, you will need to implement certain functions yourself in each 
# concrete subtype of `Polynomial`
#                                                                               
#############################################################################
#############################################################################

###################
# Polynomial type #
###################

"""
An abstract `Polynomial` type, acting as the supertype for all implementations of a polynomial.


Introduction:
The functions defined in this file, and in the folder src/basic_polynomial_operations/abstract
provide an interface for any subtype of Polynomial. Some operations (e.g. *) can 
be implemented at the abstract level assuming other functions (e.g. iterate) are implemented for
the concrete subtypes of Polynomial. This utilises Julia's ability to do multiple dispatch, where
the compiler can determine the correct function/method to call depending on the type of polynomial
it is working with.

For any function that is not implemented at the abstract level, see the PolynomialDense struct
for a particular implementation.

Note: although a function may be implemented correctly at the abstract level, without leveraging 
the details of the specific implementation the function may be very slow. Consider overriding these
functions for your concrete subtypes.


Constructors:
We assume that at least the following two constructors exist for a given concrete subtype:
    `Polynomial()` -> the empty constructor should produce the zero polynomial
    `Polynomial(vt::Vector{Term})` -> given a vector of terms `vt`, constructs a polynomial with those terms

We assume also that for any concrete subtype of Polynomial, no zero term is stored with degree higher
than the leading term (ie, we do not store 1 + 2x + x^2 + 0x^3 + 0x^7)
"""
abstract type Polynomial end

""" 
A generic error message for functionality that should be implemented by concrete subtypes.

This provides a fallback method if you forget to implement a function - ie, if Julia determines
that it cannot find the relevant function for your Polynomial subtype, it will fallback to calling
a function that throws an error.

p: 
    The polynomial for which the method cannot be applied
method: 
    The name of the unimplemented method
"""
function not_implemented_error(p::Polynomial, methodName::String)
    error("The method '$(methodName)' is not yet implemented for a polynomial of type $(typeof(p))")
end

"""
This function maintains the invariant of the Polynomial type so that there are no zero terms beyond the highest
non-zero term.
"""
function trim!(p::Polynomial)::Polynomial
    # We stop if the leading term has a non-zero coefficient or the polynomial is the zero polynomial
    while iszero(leading(p)) && !iszero(p)
        pop!(p)
    end
    return p
end

#########################
# Abstract Constructors #
#########################

"""
Construct a polynomial with a single term.
"""
function (::Type{P})(t::Term)::P where {P <: Polynomial}
    return P([t])
end

"""
Construct a polynomial of the form x^p-x.
"""
function cyclotonic_polynomial(::Type{P}, p::Int)::P where {P <: Polynomial}
    return P([Term(1,p), Term(-1,0)])
end

"""
Construct a polynomial of the form x-n.
"""
function linear_monic_polynomial(::Type{P})::P where {P <: Polynomial} 
    return P([Term(1,1), Term(-n,0)])
end

"""
Construct a polynomial of the form x.
"""
function x_poly(::Type{P})::P where {P <: Polynomial} 
    return P([Term(1,1)])
end
x_poly(p::P) where {P <: Polynomial} = x_poly(P)

"""
Creates the zero polynomial.
"""
function zero(::Type{P})::P where {P <: Polynomial} 
    return P()
end
zero(p::P) where {P <: Polynomial} = zero(P)

"""
Creates the unit polynomial.
"""
function one(::Type{P})::P where {P <: Polynomial} 
    return P(one(Term))
end
one(p::P) where {P <: Polynomial} = one(P)

"""
Generates a random polynomial.
"""
function rand(::Type{P} ; 
                degree::Int = -1, 
                terms::Int = -1, 
                max_coeff::Int = 100, 
                mean_degree::Float64 = 5.0,
                prob_term::Float64  = 0.7,
                monic = false,
                condition = (p)->true) where {P <: Polynomial}
        
    while true 
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = rand(1:max_coeff,_terms+1)
        monic && (coeffs[end] = 1)
        p = P( [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
        condition(p) && return p
    end
end

###########
# Display #
###########

"""
Show a polynomial.
"""
function show(io::IO, p::Polynomial)
    if iszero(p)
        print(io,"0")
    else
        n = length(p.terms)
        for (i,t) in enumerate(p)
            if !iszero(t)
                print(io, t, i != n ? " + " : "")
            end
        end
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.

This must be overriden by concrete subtypes.
"""
iterate(p::Polynomial, state=1) = not_implemented_error(p, "iterate")

##############################
# Queries about a polynomial #
##############################

"""
The number of (non-zero) terms of the polynomial.

This must be overriden by concrete subtypes.
"""
length(p::Polynomial) = not_implemented_error(p, "length")

"""
The leading term of the polynomial.

This must be overriden by concrete subtypes.
"""
leading(p::Polynomial) = not_implemented_error(p, "leading")

"""
Returns the coefficients of the polynomial.
"""
coeffs(p::Polynomial)::Vector{Int} = [t.coeff for t in p]

"""
The degree of the polynomial.
"""
degree(p::Polynomial)::Int = leading(p).degree 

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::Polynomial)::Int = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::Polynomial, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

################################
# Pushing and popping of terms #
################################

"""
Push a new leading term into the polynomial.
This should modify the existing polynomial p in place (ie, no new polynomials should be created).

Note: you may wish to throw an error if pushing a term of degree that is already in the 
polynomial. However, the exact implementation details are up to you.

This must be overriden by concrete subtypes.
"""
push!(p::Polynomial, t::Term) = not_implemented_error(p, "push!")

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.

Note - this should modify the existing polynomial p in place (ie, no new polynomials should be 
created).

This must be overriden by concrete subtypes.
"""
pop!(p::Polynomial)::Term = not_implemented_error(p, "pop")

"""
Check if the polynomial is zero.
"""
iszero(p::Polynomial)::Bool = iszero(leading(p)) && (degree(p) == 0)

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
-(p::P) where {P <: Polynomial} = P(map((pt)->-pt, p.terms))

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::P)::P where {P <: Polynomial} 
    der_p = P()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p, der_term)
    end
    return trim!(der_p)
end

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::Polynomial) = p ÷ content(p)

"""
A square free polynomial.
"""
square_free(p::Polynomial, prime::Int)::Polynomial = (p ÷ gcd(p,derivative(p),prime))(prime)

#############################
# Queries about polynomials #
#############################

"""
Check if two polynomials are the same.
"""
function ==(p1::Polynomial, p2::Polynomial)::Bool 
    if length(p1) != length(p2)
        return false
    end

    return all(t1 == t2 for (t1, t2) in zip(p1, p2))
end

"""
Check if a polynomial is equal to a constant `n`.
"""
function ==(p::Polynomial, n::T) where T <: Number
    degree(p) != 0 && return false
    return leading(p).coeff == n
end

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials (of the same concrete subtype).
"""
function -(p1::P, p2::P)::P where {P <: Polynomial} 
    return p1 + (-p2)
end

"""
Multiplication of polynomial and term.
"""
function *(t::Term, p1::P)::P where {P <: Polynomial} 
    return iszero(t) ? P() : P(map((pt)->t*pt, p1.terms))
end
*(p1::Polynomial, t::Term)::Polynomial = t*p1

"""
Multiplication of polynomial and an integer.
"""
*(p::Polynomial, n::Int)::Polynomial = Term(n,0)*p
*(n::Int, p::Polynomial)::Polynomial = p*n

"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::P, n::Int) where {P <: Polynomial} = (prime)->P(map((pt)->((pt ÷ n)(prime)), p.terms))

"""
Take the mod of a polynomial with an integer.
"""
function mod(f::P, p::Int)::P where {P <: Polynomial}
    f_out = P()
    for t in f
        new_t = mod(t, p)
        !iszero(new_t) && push!(f_out, new_t)
    end
    return trim!(f_out)
end

"""
Power of a polynomial mod prime.
"""
function pow_mod(p::Polynomial, n::Int, prime::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end