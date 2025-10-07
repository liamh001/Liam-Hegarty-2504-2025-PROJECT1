#############################################################################
#############################################################################
#
# A script that runs all unit tests in the project.
#                                                                               
#############################################################################
#############################################################################

using Pkg
Pkg.activate(".")

include("../poly_factorization_project.jl")

####
# Execute unit tests for integers
###
include("integers_test.jl")
println("---BEGIN INTEGER UNIT TESTS---")
test_euclid_ints()
test_ext_euclid_ints()
println("---END INTEGER UNIT TESTS---\n")

####
# Execute unit tests for data structures
###
include("heap_test.jl")
println("---BEGIN HEAP UNIT TESTS---")
test_heap()
println("---END HEAP UNIT TESTS---\n")

####
# Execute unit tests for polynomials
####
include("polynomials_test.jl")
polynomial_types = [
    PolynomialDense{Int,Int}, 
    PolynomialDense{BigInt,Int},
    PolynomialDense{ZModP{Int, 5}, Int},
    PolynomialDense{ZModP{BigInt, 5}, Int},
    PolynomialDense{ZModP{Int, 11}, Int},
    PolynomialDense{ZModP{BigInt, 11}, Int},
    PolynomialDense{ZModP{Int, 1009}, Int}]
println("---BEGIN POLYNOMIAL UNIT TESTS---\n")
for poly in polynomial_types
    println("Type of `Polynomial`: $(poly)")
    prod_test_poly(poly)
    prod_derivative_test_poly(poly)
    ext_euclid_test_poly(poly)
    division_test_poly(poly)
    println("")
end
println("---END POLYNOMIAL UNIT TESTS---\n")

####
# Execute unit tests for polynomial factorization
####
include("factorization_test.jl")
polynomial_types = [PolynomialDense{Int,Int}, PolynomialDense{BigInt,Int}]
println("---BEGIN FACTORIZATION UNIT TESTS---\n")
for poly in polynomial_types
    println("Type of `Polynomial``: $(poly)")
    factor_mod_p_test_poly(poly)
    println("")
end
println("---END FACTORISATION UNIT TESTS---")


####
# Execute unit tests for power functions
####
include("power_test.jl")
println("---BEGIN POWER FUNCTION TESTS---\n")
test_power_functions()
println("\n---END POWER FUNCTION TESTS---")
