#############################################################################
#############################################################################
#
# This file implements polynomial division for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

# # TODO - DOCSTRINGS HERE AND RENAME TO MOD_P
# """  Modular algorithm (f, g have the same concrete subtype).
# f divide by g

# f = q*g + r (mod p)

# p is a prime

# This must be overridden by concrete subtypes.
# """
# function div_mod_p(num::P, den::P) where {P <: Polynomial}
#     not_implemented_error(num, "divide")
# end

"""
The quotient from polynomial division modulo a prime. 
"""
div_mod_p(num::P, den::P, prime::Int) where {P <: Polynomial} = first(div_rem_mod_p(num, den, prime))

"""
The remainder from polynomial division modulo a prime.
"""
rem_mod_p(num::P, den::P, prime::Int) where {P <: Polynomial} = last(div_rem_mod_p(num, den, prime))
