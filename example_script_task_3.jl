using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("poly_factorization_project.jl")

using BenchmarkTools

println("=== Task 3.7: Sparse Polynomial Examples ===\n")

# Example 1: x^1000 + 1 (very sparse)
println("Example 1: x^1000 + 1")
p1_sparse = PolynomialSparse{BigInt,Int}([Term(BigInt(1), 1000), Term(BigInt(1), 0)])
q1_sparse = PolynomialSparse{BigInt,Int}([Term(BigInt(1), 1000), Term(BigInt(1), 0)])
println("Sparse: $(p1_sparse) + $(q1_sparse) = $(p1_sparse + q1_sparse)")

# Example 2: x^500 + x^250
p2_sparse = PolynomialSparse{BigInt,Int}([Term(BigInt(1), 500), Term(BigInt(1), 250)])
q2_sparse = PolynomialSparse{BigInt,Int}([Term(BigInt(2), 500), Term(BigInt(-1), 250)])
println("Sparse: $(p2_sparse) + $(q2_sparse) = $(p2_sparse + q2_sparse)")

println("\n=== Task 3.8: Performance Comparison ===\n")

# Benchmark 1: Adding sparse polynomials
println("Benchmark 1: Adding x^1000 + x^500 + 1")
sparse_p = PolynomialSparse{Int,Int}([Term(1, 1000), Term(1, 500), Term(1, 0)])
dense_p = PolynomialDense{Int,Int}([Term(1, 1000), Term(1, 500), Term(1, 0)])

println("\nSparse addition:")
@btime $sparse_p + $sparse_p
println("Dense addition:")
@btime $dense_p + $dense_p

# Benchmark 2: Multiplying sparse polynomials
println("\nBenchmark 2: Multiplying (x^100 + 1) * (x^100 + 1)")
sparse_mult = PolynomialSparse{Int,Int}([Term(1, 100), Term(1, 0)])
dense_mult = PolynomialDense{Int,Int}([Term(1, 100), Term(1, 0)])

println("\nSparse multiplication:")
@btime $sparse_mult * $sparse_mult
println("Dense multiplication:")
@btime $dense_mult * $dense_mult

# Memory usage comparison
println("\nMemory usage comparison:")
println("Sparse (x^1000 + 1): ", Base.summarysize(sparse_p), " bytes")
println("Dense (x^1000 + 1): ", Base.summarysize(dense_p), " bytes")
