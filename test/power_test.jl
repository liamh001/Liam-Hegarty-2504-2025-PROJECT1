"""
Test power functions for polynomials.
"""
function test_power_functions()
    println("=== TESTING POWER FUNCTIONS ===\n")
    Random.seed!(0)
    
    test_count = 0
    
    # Test 101 randomly generated polynomials
    for i in 1:101
        # Alternate between Dense and Sparse
        P = i % 2 == 1 ? PolynomialDense{Int,Int} : PolynomialSparse{Int,Int}
        
        p = rand(P; degree=3, max_coeff=10)
        n = rand(2:8)
        
        # Test ^ using repeated squaring vs naive multiplication
        result = p^n
        expected = p
        for j in 2:n
            expected *= p
        end
        @assert result == expected
        
        # Test pow_mod
        prime = rand([5, 7, 11, 13])
        result_mod = pow_mod(p, n, prime)
        expected_mod = mod(expected, prime)
        @assert mod(result_mod - expected_mod, prime) == 0
        
        test_count += 1
        
        if i % 20 == 0
            print(".")
        end
    end
    
    println("\n✓ All 101 random polynomial tests passed")
    println("✓ Both ^ and pow_mod working correctly with repeated squaring")
    println("\n=== ALL POWER TESTS PASSED ===")
end
