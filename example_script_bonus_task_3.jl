#############################################################################
#############################################################################
#
# Example Script for Bonus Task 3: Pseudo-Division in Z[x]
# 
# This script demonstrates pseudo-division over Z[x] using BigInt
# to handle coefficient swell.
#
#############################################################################

using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("poly_factorization_project.jl")

#############################################################################
# PSEUDO-DIVISION IMPLEMENTATION
#############################################################################

"""
Pseudo-division for polynomials over Z[x].
"""
function pseudo_div_rem(f::P, g::P)::Tuple{P, P} where {C <: Integer, D, P <: Polynomial{C, D}}
    iszero(g) && error("Cannot divide by zero polynomial")
    
    if degree(f) < degree(g)
        lc_g = leading(g).coeff
        return (zero(P), lc_g * f)
    end
    
    q = zero(P)
    r = deepcopy(f)
    c = leading(g).coeff
    deg_g = degree(g)
    
    while degree(r) >= deg_g && !iszero(r)
        lc_r = leading(r).coeff
        deg_r = degree(r)
        
        # Multiply both q and r by c
        q = c * q
        r = c * r
        
        # Create and subtract the term
        s = P(Term(lc_r, deg_r - deg_g))
        q = q + s
        r = r - s * g
        
        r = trim!(r)
    end
    
    return (trim!(q), trim!(r))
end

function pseudo_quo(f::P, g::P)::P where {P <: Polynomial}
    return first(pseudo_div_rem(f, g))
end

function pseudo_rem(f::P, g::P)::P where {P <: Polynomial}
    return last(pseudo_div_rem(f, g))
end

function pseudo_gcd(f::P, g::P) where {P <: Polynomial}
    while !(iszero(f) || iszero(g))
        f, g = g, pseudo_rem(f, g)
    end
    iszero(g) && return prim_part(f)
    return prim_part(g)
end

function square_free_pseudo(f::P) where {C <: Integer, D, P <: Polynomial{C, D}}
    f = prim_part(f)
    df = derivative(f)
    g = pseudo_gcd(f, df)
    
    if iszero(g) || degree(g) == 0
        return f
    end
    
    sq_fr_f = pseudo_quo(f, g)
    return prim_part(sq_fr_f)
end

#############################################################################
# DEMONSTRATIONS
#############################################################################

println("="^80)
println("BONUS TASK 3: PSEUDO-DIVISION IN Z[x]")
println("="^80)

println("\n" * "="^80)
println("Part 1: Pseudo-Division Property")
println("="^80)

# Use BigInt to avoid overflow
x = x_poly(PolynomialDense{BigInt, Int})

println("\nThe pseudo-division property states:")
println("Given f, g in Z[x], pseudo_div_rem(f, g) returns (q, r) such that:")
println("    lc(g)^delta * f = q * g + r")
println("where delta is at most deg(f) - deg(g) + 1 and deg(r) < deg(g)\n")

# Example 1
println("Example 1:")
println("-"^50)
f1 = 3*x^3 + 2*x^2 + x + 1
g1 = 2*x^2 + x
println("f = ", f1)
println("g = ", g1)

q1, r1 = pseudo_div_rem(f1, g1)
println("Quotient q = ", q1)
println("Remainder r = ", r1)

# Verify property
println("\nVerification: q * g + r should equal some power of lc(g) times f")
lhs1 = q1 * g1 + r1
println("q * g + r = ", lhs1)
println("Degree of r < degree of g? ", degree(r1) < degree(g1))

# Example 2
println("\n\nExample 2:")
println("-"^50)
f2 = 5*x^3 + 4*x^2 + 3*x + 2
g2 = 2*x + 1
println("f = ", f2)
println("g = ", g2)

q2, r2 = pseudo_div_rem(f2, g2)
println("Quotient q = ", q2)
println("Remainder r = ", r2)
println("Degree of r < degree of g? ", degree(r2) < degree(g2))

# Example 3
println("\n\nExample 3:")
println("-"^50)
f3 = x^4 + x^3 + x^2 + x + 1
g3 = x^2 + 1
println("f = ", f3)
println("g = ", g3)

q3, r3 = pseudo_div_rem(f3, g3)
println("Quotient q = ", q3)
println("Remainder r = ", r3)
println("Degree of r < degree of g? ", degree(r3) < degree(g3))

#############################################################################
# Part 2: Pseudo-GCD
#############################################################################

println("\n" * "="^80)
println("Part 2: Pseudo-GCD Property")
println("="^80)

println("\nThe pseudo-gcd divides both polynomials.")

# Example 1
println("\nExample 1:")
println("-"^50)
p1 = (x + 1) * (x + 2)
p2 = (x + 1) * (x + 3)
println("p1 = ", p1)
println("p2 = ", p2)

gcd1 = pseudo_gcd(p1, p2)
println("pseudo_gcd = ", gcd1)

r1a = pseudo_rem(p1, gcd1)
r1b = pseudo_rem(p2, gcd1)
println("pseudo_rem(p1, gcd) = ", r1a)
println("pseudo_rem(p2, gcd) = ", r1b)
println("Both zero? ", iszero(r1a) && iszero(r1b))

# Example 2
println("\n\nExample 2:")
println("-"^50)
p3 = (x^2 + 1) * (x + 2)
p4 = (x^2 + 1) * (x + 3)
println("p3 = ", p3)
println("p4 = ", p4)

gcd2 = pseudo_gcd(p3, p4)
println("pseudo_gcd = ", gcd2)

r2a = pseudo_rem(p3, gcd2)
r2b = pseudo_rem(p4, gcd2)
println("Both zero? ", iszero(r2a) && iszero(r2b))

#############################################################################
# Part 3: Square-Free Factorization
#############################################################################

println("\n" * "="^80)
println("Part 3: Square-Free Factorization")
println("="^80)

primes = [5, 7, 11]

function make_monic_mod_p(p::P, prime::Int) where {P <: Polynomial}
    p_mod = mod(p, prime)
    lc = leading(p_mod).coeff
    if lc == 0
        return p_mod
    end
    inv_lc = int_inverse_mod(Int(lc), prime)
    return mod(p_mod * inv_lc, prime)
end

# Example 1
println("\nExample 1: (x + 1)^2 * (x + 2)")
println("-"^50)
sf1 = (x + 1)^2 * (x + 2)
println("Original: ", sf1)

sf1_sqfree = square_free_pseudo(sf1)
println("Square-free: ", sf1_sqfree)

for prime in primes
    sf1_mod = make_monic_mod_p(sf1_sqfree, prime)
    sf1_direct = square_free_mod_p(sf1, prime)
    sf1_direct = make_monic_mod_p(sf1_direct, prime)
    println("Prime ", prime, ": Equal? ", sf1_mod == sf1_direct)
end

# Example 2
println("\n\nExample 2: (x + 2)^3 * (x + 3)")
println("-"^50)
sf2 = (x + 2)^3 * (x + 3)
println("Original: ", sf2)

sf2_sqfree = square_free_pseudo(sf2)
println("Square-free: ", sf2_sqfree)

for prime in primes
    sf2_mod = make_monic_mod_p(sf2_sqfree, prime)
    sf2_direct = square_free_mod_p(sf2, prime)
    sf2_direct = make_monic_mod_p(sf2_direct, prime)
    println("Prime ", prime, ": Equal? ", sf2_mod == sf2_direct)
end

# Example 3
println("\n\nExample 3: (x + 1)^2 * (x + 3)^2")
println("-"^50)
sf3 = (x + 1)^2 * (x + 3)^2
println("Original: ", sf3)

sf3_sqfree = square_free_pseudo(sf3)
println("Square-free: ", sf3_sqfree)

for prime in primes
    sf3_mod = make_monic_mod_p(sf3_sqfree, prime)
    sf3_direct = square_free_mod_p(sf3, prime)
    sf3_direct = make_monic_mod_p(sf3_direct, prime)
    println("Prime ", prime, ": Equal? ", sf3_mod == sf3_direct)
end

# Example 4
println("\n\nExample 4: (x^2 + 1)^2 * (x + 1)")
println("-"^50)
sf4 = (x^2 + 1)^2 * (x + 1)
println("Original: ", sf4)

sf4_sqfree = square_free_pseudo(sf4)
println("Square-free: ", sf4_sqfree)

for prime in primes
    sf4_mod = make_monic_mod_p(sf4_sqfree, prime)
    sf4_direct = square_free_mod_p(sf4, prime)
    sf4_direct = make_monic_mod_p(sf4_direct, prime)
    println("Prime ", prime, ": Equal? ", sf4_mod == sf4_direct)
end

# Example 5
println("\n\nExample 5: (x + 1)^4 * (x + 2)^2")
println("-"^50)
sf5 = (x + 1)^4 * (x + 2)^2
println("Original: ", sf5)

sf5_sqfree = square_free_pseudo(sf5)
println("Square-free: ", sf5_sqfree)

for prime in primes
    sf5_mod = make_monic_mod_p(sf5_sqfree, prime)
    sf5_direct = square_free_mod_p(sf5, prime)
    sf5_direct = make_monic_mod_p(sf5_direct, prime)
    println("Prime ", prime, ": Equal? ", sf5_mod == sf5_direct)
end

println("\n" * "="^80)
println("BONUS TASK 3 COMPLETE")
println("="^80)

println("\n\n" * "="^80)
println("TECHNICAL EXPLANATION")
println("="^80)

println("""
PSEUDO-DIVISION ALGORITHM

In Z[x], we cannot always perform exact division because Z is not a field.
For example, dividing 3x + 2 by 2x + 1 would require fractional coefficients.

Pseudo-division multiplies by powers of the divisor's leading coefficient
to ensure all operations stay in Z[x].

Given f, g in Z[x], compute q, r such that:
    c^k * f = q * g + r  for some power k
where c = lc(g) and deg(r) < deg(g).

The major problem is COEFFICIENT SWELL - coefficients grow exponentially.
This is why we use BigInt in this implementation.

In Q[x] (rational coefficients), we have exact division without coefficient
swell. The trade-off is storing fractions, but this is more efficient than
the exponential growth in Z[x].

Pseudo-division works but is impractical for large computations due to
coefficient growth. Modern systems use Q[x] or modular methods instead.
""")#
