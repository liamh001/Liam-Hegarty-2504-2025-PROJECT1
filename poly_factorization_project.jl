#############################################################################
#############################################################################
#
# This is the main project file for polynomial factorization
#                                                                               
#############################################################################
#############################################################################

using Distributions, StatsBase, Random

import Base: %
import Base: push!, pop!, iszero, show, isless, map, map!, iterate, length, last
import Base: +, -, *, mod, %, ÷, ==, ^, rand, rem, zero, one

include("src/general_alg.jl")
include("src/term.jl")

# Polynomial
include("src/polynomial_definitions/polynomial.jl")
    include("src/basic_polynomial_operations/abstract/polynomial_addition.jl")
    include("src/basic_polynomial_operations/abstract/polynomial_multiplication.jl")
    include("src/basic_polynomial_operations/abstract/polynomial_division.jl")
    include("src/basic_polynomial_operations/abstract/polynomial_gcd.jl")

# Dense
include("src/polynomial_definitions/polynomial_dense.jl")
    include("src/basic_polynomial_operations/dense/polynomial_addition.jl")
    include("src/basic_polynomial_operations/dense/polynomial_multiplication.jl")
    include("src/basic_polynomial_operations/dense/polynomial_division.jl")
    include("src/basic_polynomial_operations/dense/polynomial_gcd.jl")

include("src/polynomial_factorization/factor.jl")

nothing