#############################################################################
#############################################################################
#
# This file contains units tests for polynomial factorization
#                                                                               
#############################################################################
#############################################################################


"""
Test factorization of polynomials.
"""
function factor_test_poly(::Type{P};
    N::Int = 5, seed::Int = 0, primes::Vector{Int} = [5,11,13]
    ) where {P <: Polynomial}
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(P; max_coeff = 3, mean_degree = 2.0)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end

    println("\nfactor_test_poly for $(P) - PASSED")
end

