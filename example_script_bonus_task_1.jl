#############################################################################
#
# Example Script for Bonus Task 1: CRT Multiplication
# 
# This reuses the integer CRT and polynomial CRT from Task 7.
#
#############################################################################

using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("poly_factorization_project.jl")
using BenchmarkTools

println("="^80)
println("BONUS TASK 1: CRT MULTIPLICATION")
println("="^80)

#############################################################################
# Bonus Task 1.1: Implementation
#############################################################################

"""
Multiply two polynomials using Chinese Remainder Theorem.

Algorithm:
    Input: f, g in Z[x]
    Output: f * g in Z[x]
    
    1. Compute coefficient bound: B = 2 * height(f) * height(g) * min(deg(f)+1, deg(g)+1)
       This bounds the maximum coefficient in the product f*g
    2. Select primes p1, p2, ..., pk such that p1*p2*...*pk > B
    3. For each prime pi: compute ci = (f * g) mod pi
    4. Use CRT to reconstruct c from c1, c2, ..., ck
    5. Apply symmetric mod to get final result in Z[x]
"""
function crt_multiplication(f::P, g::P) where {C, D, P <: Polynomial{C, D}}
    # Handle trivial cases
    iszero(f) && return zero(P)
    iszero(g) && return zero(P)
    
    # Compute coefficient bound for the product
    # The largest coefficient in f*g is bounded by:
    # 2 * max(|f_i|) * max(|g_j|) * min(deg(f)+1, deg(g)+1)
    height_f = poly_height(f)
    height_g = poly_height(g)
    n = degree(f)
    m = degree(g)
    
    B = 2 * height_f * height_g * min(n + 1, m + 1)
    
    # Select primes until their product exceeds B
    primes = Int[]
    M = BigInt(1)
    p = 2
    
    while M <= B
        p = next_prime(p)
        # Skip primes that divide leading coefficients (bad primes)
        if mod(leading(f).coeff, p) != 0 && mod(leading(g).coeff, p) != 0
            push!(primes, p)
            M *= p
        end
    end
    
    # Compute product modulo each prime
    products_mod_p = P[]
    
    for prime in primes
        f_mod = mod(f, prime)
        g_mod = mod(g, prime)
        prod_mod = mod(f_mod * g_mod, prime)
        push!(products_mod_p, prod_mod)
    end
    
    # Reconstruct using CRT (int_crt already applies symmetric mod)
    result = poly_crt(products_mod_p, primes)
    
    # Just trim to remove any zero leading terms
    return trim!(result)
end

#############################################################################
# Bonus Task 1.2: Correctness Tests (10+ non-trivial examples)
#############################################################################

println("\n" * "="^80)
println("BONUS TASK 1.2: CORRECTNESS TESTS")
println("="^80)

x = x_poly(PolynomialDense{BigInt, Int})

# Test 1: Using binomial theorem (x + 1)^n
println("\nTest 1: (x + 1)^50 * (x + 1)^51 = (x + 1)^101")
p1 = (x + BigInt(1))^50
p2 = (x + BigInt(1))^51
expected1 = (x + BigInt(1))^101
result1 = crt_multiplication(p1, p2)
println("  Degree: ", degree(result1), " (expected ", degree(expected1), ")")
println("  Match: ", result1 == expected1 ? "PASS" : "FAIL")

# Test 2: (2x + 3)^n using binomial theorem
println("\nTest 2: (2x + 3)^60 * (2x + 3)^45 = (2x + 3)^105")
p3 = (2*x + BigInt(3))^60
p4 = (2*x + BigInt(3))^45
expected2 = (2*x + BigInt(3))^105
result2 = crt_multiplication(p3, p4)
println("  Degree: ", degree(result2), " (expected ", degree(expected2), ")")
println("  Match: ", result2 == expected2 ? "PASS" : "FAIL")

# Test 3: Difference of squares pattern
println("\nTest 3: (x^50 + 1) * (x^50 - 1) = x^100 - 1")
p5 = x^50 + BigInt(1)
p6 = x^50 + BigInt(-1)  # Change here
expected3 = x^100 + BigInt(-1)  # Change here
result3 = crt_multiplication(p5, p6)
println("  Degree: ", degree(result3), " (expected ", degree(expected3), ")")
println("  Match: ", result3 == expected3 ? "PASS" : "FAIL")


# Test 4: Large powers
println("\nTest 4: (x + 5)^75 * (x + 5)^75 = (x + 5)^150")
p7 = (x + BigInt(5))^75
expected4 = (x + BigInt(5))^150
result4 = crt_multiplication(p7, p7)
println("  Degree: ", degree(result4), " (expected ", degree(expected4), ")")
println("  Match: ", result4 == expected4 ? "PASS" : "FAIL")

# Test 5: Sparse polynomial
println("\nTest 5: (x^100 + x^50 + 1) * (x^100 - x^50 + 1)")
p8 = x^100 + x^50 + BigInt(1)
p9 = x^100 - x^50 + BigInt(1)
expected5 = p8 * p9
result5 = crt_multiplication(p8, p9)
println("  Degree: ", degree(result5), " (expected ", degree(expected5), ")")
println("  Match: ", result5 == expected5 ? "PASS" : "FAIL")

# Test 6: Large coefficients
println("\nTest 6: (1000000x + 999999)^55 * (1000000x + 999999)^55")
p10 = (BigInt(1000000)*x + BigInt(999999))^55
expected6 = (BigInt(1000000)*x + BigInt(999999))^110
result6 = crt_multiplication(p10, p10)
println("  Degree: ", degree(result6), " (expected ", degree(expected6), ")")
println("  Match: ", result6 == expected6 ? "PASS" : "FAIL")

# Test 7: Mixed terms
println("\nTest 7: (3x + 7)^80 * (3x + 7)^30 = (3x + 7)^110")
p11 = (3*x + BigInt(7))^80
p12 = (3*x + BigInt(7))^30
expected7 = (3*x + BigInt(7))^110
result7 = crt_multiplication(p11, p12)
println("  Degree: ", degree(result7), " (expected ", degree(expected7), ")")
println("  Match: ", result7 == expected7 ? "PASS" : "FAIL")

# Test 8: Quadratic base
println("\nTest 8: (x^2 + x + 1)^60 * (x^2 + x + 1)^60")
p13 = (x^2 + x + BigInt(1))^60
expected8 = (x^2 + x + BigInt(1))^120
result8 = crt_multiplication(p13, p13)
println("  Degree: ", degree(result8), " (expected ", degree(expected8), ")")
println("  Match: ", result8 == expected8 ? "PASS" : "FAIL")

# Test 9: Negative coefficients
println("\nTest 9: (x - 10)^90 * (x - 10)^20 = (x - 10)^110")
p14 = (x + BigInt(-10))^90
p15 = (x + BigInt(-10))^20
expected9 = (x + BigInt(-10))^110
result9 = crt_multiplication(p14, p15)
println("  Degree: ", degree(result9), " (expected ", degree(expected9), ")")
println("  Match: ", result9 == expected9 ? "PASS" : "FAIL")

# Test 10: Larger degree base
println("\nTest 10: (x^3 + 2x^2 + 3x + 4)^40 * (x^3 + 2x^2 + 3x + 4)^40")
p16 = (x^3 + 2*x^2 + 3*x + BigInt(4))^40
expected10 = (x^3 + 2*x^2 + 3*x + BigInt(4))^80
result10 = crt_multiplication(p16, p16)
println("  Degree: ", degree(result10), " (expected ", degree(expected10), ")")
println("  Match: ", result10 == expected10 ? "PASS" : "FAIL")

# Test 11: Very sparse polynomial
println("\nTest 11: (x^150 + 1) * (x^150 + 1) = x^300 + 2x^150 + 1")
p17 = x^150 + BigInt(1)
expected11 = x^300 + 2*x^150 + BigInt(1)
result11 = crt_multiplication(p17, p17)
println("  Degree: ", degree(result11), " (expected ", degree(expected11), ")")
println("  Match: ", result11 == expected11 ? "PASS" : "FAIL")

# Test 12: Different bases multiplied
println("\nTest 12: (x + 2)^70 * (x + 3)^70")
p18 = (x + BigInt(2))^70
p19 = (x + BigInt(3))^70
expected12 = p18 * p19
result12 = crt_multiplication(p18, p19)
println("  Degree: ", degree(result12), " (expected ", degree(expected12), ")")
println("  Match: ", result12 == expected12 ? "PASS" : "FAIL")

#############################################################################
# Bonus Task 1.3: Benchmarking - PolynomialDense
#############################################################################

println("\n" * "="^80)
println("BONUS TASK 1.3: BENCHMARKING - PolynomialDense with BigInt")
println("="^80)

# Benchmark 1: Moderate degree
println("\nBenchmark 1: (x + 1)^60 * (x + 1)^60")
b_p1 = (x + BigInt(1))^60
println("Direct multiplication:")
@btime $b_p1 * $b_p1
println("CRT multiplication:")
@btime crt_multiplication($b_p1, $b_p1)

# Benchmark 2: Higher degree
println("\nBenchmark 2: (x + 2)^80 * (x + 3)^80")
b_p2 = (x + BigInt(2))^80
b_p3 = (x + BigInt(3))^80
println("Direct multiplication:")
@btime $b_p2 * $b_p3
println("CRT multiplication:")
@btime crt_multiplication($b_p2, $b_p3)

# Benchmark 3: Larger coefficients
println("\nBenchmark 3: (10000x + 9999)^50 * (10000x + 9999)^50")
b_p4 = (BigInt(10000)*x + BigInt(9999))^50
println("Direct multiplication:")
@btime $b_p4 * $b_p4
println("CRT multiplication:")
@btime crt_multiplication($b_p4, $b_p4)

#############################################################################
# Bonus Task 1.3: Benchmarking - PolynomialSparse
#############################################################################

println("\n" * "="^80)
println("BONUS TASK 1.3: BENCHMARKING - PolynomialSparse with BigInt")
println("="^80)

x_sparse = x_poly(PolynomialSparse{BigInt, Int})

# Benchmark 4: Sparse polynomials
println("\nBenchmark 4: (x^100 + 1) * (x^100 + 1)")
b_p5 = x_sparse^100 + BigInt(1)
println("Direct multiplication:")
@btime $b_p5 * $b_p5
println("CRT multiplication:")
@btime crt_multiplication($b_p5, $b_p5)

# Benchmark 5: More terms in sparse polynomial
println("\nBenchmark 5: (x^150 + x^75 + 1) * (x^150 - x^75 + 1)")
b_p6 = x_sparse^150 + x_sparse^75 + BigInt(1)
b_p7 = x_sparse^150 - x_sparse^75 + BigInt(1)
println("Direct multiplication:")
@btime $b_p6 * $b_p7
println("CRT multiplication:")
@btime crt_multiplication($b_p6, $b_p7)

# Benchmark 6: Very sparse, high degree
println("\nBenchmark 6: (x^200 + x^100 + x^50 + 1)^2")
b_p8 = x_sparse^200 + x_sparse^100 + x_sparse^50 + BigInt(1)
println("Direct multiplication:")
@btime $b_p8 * $b_p8
println("CRT multiplication:")
@btime crt_multiplication($b_p8, $b_p8)
