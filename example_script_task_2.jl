using Pkg
Pkg.activate(".")  # Activate the project in current directory
Pkg.instantiate()
include("poly_factorization_project.jl")

# Example 1: typemax(Int) + 1
p1 = PolynomialDense{Int,Int}([Term(typemax(Int), 0), Term(0, 1)])
q1 = PolynomialDense{Int,Int}([Term(1, 0)])
println("Example 1: $(typemax(Int)) + 1 = $((p1 + q1).terms[1].coeff)")

# Example 2: Large number addition
p2 = PolynomialDense{Int,Int}([Term(9223372036854775800, 0)])
q2 = PolynomialDense{Int,Int}([Term(100, 0)])
println("Example 2: 9223372036854775800 + 100 = $((p2 + q2).terms[1].coeff)")

# Example 3: 5×10^18 + 6×10^18
p3 = PolynomialDense{Int,Int}([Term(5000000000000000000, 0), Term(0, 1), Term(1, 2)])
q3 = PolynomialDense{Int,Int}([Term(6000000000000000000, 0), Term(0, 1), Term(2, 2)])
result3 = p3 + q3
println("Example 3: 5×10^18 + 6×10^18 = $(result3.terms[1].coeff)")

# Example 4: typemin(Int) - 1
p4 = PolynomialDense{Int,Int}([Term(typemin(Int), 0)])
q4 = PolynomialDense{Int,Int}([Term(-1, 0)])
println("Example 4: $(typemin(Int)) - 1 = $((p4 + q4).terms[1].coeff)")

# Example 5: Both coefficients overflow
p5 = PolynomialDense{Int,Int}([Term(4611686018427387904, 0), Term(4611686018427387904, 1)])
q5 = PolynomialDense{Int,Int}([Term(4611686018427387904, 0), Term(4611686018427387904, 1)])
result5 = p5 + q5
println("Example 5: Both coefficients overflow: $(result5)")

#after parametrising poly

# Example 1 with BigInt
p1_big = PolynomialDense{BigInt,Int}([Term(BigInt(typemax(Int)), 0), Term(BigInt(0), 1)])
q1_big = PolynomialDense{BigInt,Int}([Term(BigInt(1), 0)])
println("Example 1 (BigInt): $(typemax(Int)) + 1 = $((p1_big + q1_big).terms[1].coeff)")

# Example 2 with BigInt
p2_big = PolynomialDense{BigInt,Int}([Term(BigInt(9223372036854775800), 0)])
q2_big = PolynomialDense{BigInt,Int}([Term(BigInt(100), 0)])
println("Example 2 (BigInt): 9223372036854775800 + 100 = $((p2_big + q2_big).terms[1].coeff)")

# Example 3 with BigInt
p3_big = PolynomialDense{BigInt,Int}([Term(BigInt(5000000000000000000), 0), Term(BigInt(0), 1), Term(BigInt(1), 2)])
q3_big = PolynomialDense{BigInt,Int}([Term(BigInt(6000000000000000000), 0), Term(BigInt(0), 1), Term(BigInt(2), 2)])
result3_big = p3_big + q3_big
println("Example 3 (BigInt): 5×10^18 + 6×10^18 = $(result3_big.terms[1].coeff)")

# Example 4 with BigInt
p4_big = PolynomialDense{BigInt,Int}([Term(BigInt(typemin(Int)), 0)])
q4_big = PolynomialDense{BigInt,Int}([Term(BigInt(-1), 0)])
println("Example 4 (BigInt): $(typemin(Int)) - 1 = $((p4_big + q4_big).terms[1].coeff)")

# Example 5 with BigInt
p5_big = PolynomialDense{BigInt,Int}([Term(BigInt(4611686018427387904), 0), Term(BigInt(4611686018427387904), 1)])
q5_big = PolynomialDense{BigInt,Int}([Term(BigInt(4611686018427387904), 0), Term(BigInt(4611686018427387904), 1)])
result5_big = p5_big + q5_big
println("Example 5 (BigInt): Both coefficients correct: $(result5_big)")