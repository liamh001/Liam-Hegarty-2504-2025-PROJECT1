using Pkg
Pkg.activate(".")
include("poly_factorization_project.jl")

println("=== Task 5: ZModP Examples ===\n")

# Example with p = 5
println("Working in Z_5:")
a = ZModP(3, Val(5))
b = ZModP(4, Val(5))

println("a = $a, b = $b")
println("a + b = $(a + b) (3 + 4 = 7 ≡ 2 mod 5)")
println("a * b = $(a * b) (3 * 4 = 12 ≡ 2 mod 5)")
println("a - b = $(a - b) (3 - 4 = -1 ≡ 4 mod 5)")
println("a ÷ b = $(a ÷ b) (3 * 4^(-1) = 3 * 4 = 12 ≡ 2 mod 5)")
println("a^3 = $(a^3) (3^3 = 27 ≡ 2 mod 5)")
println("-a = $(-a) (-3 ≡ 2 mod 5)")
println("inv(a) = $(inv(a)) (3^(-1) ≡ 2 mod 5)")

# Example with p = 7
println("\nWorking in Z_7:")
c = ZModP(5, Val(7))
println("c = $c")
println("c + 3 = $(c + 3)")
println("2 * c = $(2 * c)")
println("c^6 = $(c^6) (Fermat's little theorem: 5^6 ≡ 1 mod 7)")

# Example with polynomials
println("\nPolynomial example:")
x = x_poly(PolynomialDense{ZModP{Int, 7}, Int})
f = x^2 + 5*x + 1
println("f = $f in Z_7[x]")
println("f(2) = $(evaluate(f, ZModP(2, Val(7))))")

# Demonstrate error handling
println("\nError handling:")
try
    bad = ZModP(5, Val(6))  # 6 is not prime
catch e
    println("Expected error: ", e)
end
