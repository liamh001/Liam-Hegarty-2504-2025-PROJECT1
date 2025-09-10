#############################################################################
#############################################################################
#
# This file implements polynomial division for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

# TODO - DOCSTRINGS HERE AND RENAME TO MOD_P
"""  Modular algorithm (f, g have the same concrete subtype).
f divide by g

f = q*g + r (mod p)

p is a prime

This must be overridden by concrete subtypes.
"""
function divide(num::P, den::P) where {P <: Polynomial}
    not_implemented_error(num, "divide")
end

"""
The quotient from polynomial division. Returns a function of an integer.
"""
÷(num::P, den::P) where {P <: Polynomial}  = (p::Int) -> first(divide(num,den)(p))

"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::P, den::P) where {P <: Polynomial} = (p::Int) -> last(divide(num,den)(p))
