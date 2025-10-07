#############################################################################
#
# Example Script for Task 7: Chinese Remainder Theorem
# 
# This script demonstrates:
# 1. Integer CRT with 5+ remainders/moduli (3+ examples)
# 2. Polynomial CRT 
# 3. Polynomial multiplication using CRT
#
#############################################################################

# Activate the project environment and install dependencies
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# Load the polynomial project
include("poly_factorization_project.jl")

println("="^70)
println("Task 7: Chinese Remainder Theorem Examples")
println("="^70)

#############################################################################
# Integer CRT Examples (demonstrating correctness)
#############################################################################

println("\n" * "="^70)
println("Part 1: Integer CRT Examples (5+ remainders/moduli)")
println("="^70)

# Example 1: From the document (extended to 5 moduli)
println("\nExample 1:")
println("-"^50)
rems1 = [2, 1, 5, 3, 4]
mods1 = [3, 5, 7, 11, 13]
println("Remainders: ", rems1)
println("Moduli:     ", mods1)
result1 = int_crt(rems1, mods1)
println("Result:     ", result1)
println("\nVerification:")
for i in 1:length(rems1)
    computed = mod(result1, mods1[i])
    expected = rems1[i]
    match = computed == expected ? "pass" : "FAIL"
    println("  $result1 mod $(mods1[i]) = $computed (expected $(expected)) $match")
end

# Example 2: Larger numbers
println("\nExample 2:")
println("-"^50)
rems2 = [15, 22, 8, 31, 17, 44]
mods2 = [17, 23, 29, 37, 41, 47]
println("Remainders: ", rems2)
println("Moduli:     ", mods2)
result2 = int_crt(rems2, mods2)
println("Result:     ", result2)
println("\nVerification:")
for i in 1:length(rems2)
    computed = mod(result2, mods2[i])
    expected = rems2[i]
    match = computed == expected ? "pass" : "FAIL"
    println("  $result2 mod $(mods2[i]) = $computed (expected $(expected)) $match")
end

# Example 3: Starting from a known number
println("\nExample 3:")
println("-"^50)
original_number = 123456789
mods3 = [7, 11, 13, 17, 19, 23, 29]
rems3 = [mod(original_number, m) for m in mods3]
println("Original number: ", original_number)
println("Moduli:          ", mods3)
println("Remainders:      ", rems3)
result3 = int_crt(rems3, mods3)
println("Reconstructed:   ", result3)
M_total = prod(BigInt.(mods3))
# Result should equal original_number mod M_total
expected = mod(BigInt(original_number), M_total)
if expected > M_total / 2
    expected -= M_total
end
println("Expected:        ", expected)
println("Match:           ", result3 == expected ? "pass" : "FAIL")

#############################################################################
# Polynomial CRT Examples
#############################################################################

println("\n" * "="^70)
println("Part 2: Polynomial CRT Examples")
println("="^70)

x = x_poly(PolynomialDense{Int, Int})

# Example from the document: (3x - 4)(6x + 5) = 18x^2 - 9x - 20
println("\nExample: Reconstructing (3x - 4)(6x + 5) using CRT")
println("-"^50)
a = 3*x + (-4)
b = 6*x + 5
c_expected = a * b
println("a = ", a)
println("b = ", b)
println("Expected product c = ", c_expected)

# Factor modulo different primes (from document table)
p1 = 5
p2 = 7
p3 = 3

a_p1 = mod(a, p1)
b_p1 = mod(b, p1)
c_p1 = mod(a_p1 * b_p1, p1)
println("\nmod $p1: a = $a_p1, b = $b_p1, c = $c_p1")

a_p2 = mod(a, p2)
b_p2 = mod(b, p2)
c_p2 = mod(a_p2 * b_p2, p2)
println("mod $p2: a = $a_p2, b = $b_p2, c = $c_p2")

a_p3 = mod(a, p3)
b_p3 = mod(b, p3)
c_p3 = mod(a_p3 * b_p3, p3)
println("mod $p3: a = $a_p3, b = $b_p3, c = $c_p3")

# Reconstruct using CRT
c_reconstructed = poly_crt([c_p1, c_p2, c_p3], [p1, p2, p3])
println("\nReconstructed c = ", c_reconstructed)
println("Expected c      = ", c_expected)
println("Match:            ", c_reconstructed == c_expected ? "pass" : "FAIL")

