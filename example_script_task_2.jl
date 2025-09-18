using Pkg
Pkg.activate(".")  # Activate the project in current directory
Pkg.instantiate()

include("poly_factorization_project.jl")



p1 = PolynomialDense{Int,Int}([Term(typemax(Int), 0), Term(0, 1)])
q1 = PolynomialDense{Int,Int}([Term(1, 0)])
println("Example 1: $(typemax(Int)) + 1 = $((p1 + q1).terms[1].coeff)")


p2 = PolynomialDense([9223372036854775800])
q2 = PolynomialDense([100])
println("Example 2: 9223372036854775800 + 100 = $((p2 + q2).terms[1].coeff)")

p3 = PolynomialDense([5000000000000000000, 0, 1])
q3 = PolynomialDense([6000000000000000000, 0, 2])
result3 = p3 + q3
println("Example 3: 5×10^18 + 6×10^18 = $(result3.terms[1].coeff)")

p4 = PolynomialDense([typemin(Int)])
q4 = PolynomialDense([-1])
println("Example 4: $(typemin(Int)) - 1 = $((p4 + q4).terms[1].coeff)")

p5 = PolynomialDense([4611686018427387904, 4611686018427387904])
q5 = PolynomialDense([4611686018427387904, 4611686018427387904])
result5 = p5 + q5
println("Example 5: Both coefficients overflow: $(result5)")
