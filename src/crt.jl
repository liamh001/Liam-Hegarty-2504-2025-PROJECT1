#############################################################################
#############################################################################
#
# This file implements Chinese Remainder Theorem for Task 7
# Task 7: Factoring with Chinese Remainder Theorem over Z[x]
#
#############################################################################
#############################################################################

"""
Chinese Remainder Theorem for integers.


Returns: x mod (m1 * m2 * ... * mk)

Implementation follows the constructive proof:
    u = v1 + v2*m1 + v3*m1*m2 + ... + vk*m1*m2*...*m(k-1)
where:
    v1 <- u1
    v2 <- (u2 - v1) * (m1)^(-1) mod m2
    v3 <- (u3 - v1 - v2*m1) * (m1*m2)^(-1) mod m3
    ...
"""
function int_crt(rems::Vector{T}, moduli::Vector{T}) where {T <: Union{Int, BigInt}}
    @assert length(rems) == length(moduli) > 0
    
    k = length(rems)
    
    # Use BigInt for intermediate calculations to avoid overflow
    v = Vector{BigInt}(undef, k)
    
    # v1 <- u1
    v[1] = BigInt(rems[1])
    
    # Compute remaining v coefficients
    for i in 2:k
        # Compute the sum: v1 + v2*m1 + v3*m1*m2 + ... + v(i-1)*m1*m2*...*m(i-2)
        sum_term = v[1]
        M_prod = BigInt(1)
        for j in 2:(i-1)
            M_prod *= BigInt(moduli[j-1])
            sum_term += v[j] * M_prod
        end
        
        # vi <- (ui - sum_term) * (m1*m2*...*m(i-1))^(-1) mod mi
        M_prod *= BigInt(moduli[i-1])  # Now M_prod = m1*m2*...*m(i-1)
        m_i = BigInt(moduli[i])
        
        v[i] = mod((BigInt(rems[i]) - sum_term) * int_inverse_mod(M_prod, m_i), m_i)
    end
    
    # Construct result: u = v1 + v2*m1 + v3*m1*m2 + ...
    result = v[1]
    M_prod = BigInt(1)
    for i in 2:k
        M_prod *= BigInt(moduli[i-1])
        result += v[i] * M_prod
    end
    
    # Apply symmetric mod
    M_total = BigInt(1)
    for m in moduli
        M_total *= m
    end
    
    result = mod(result, M_total)
    if result > M_total / 2
        result -= M_total
    end
    
    # Return as BigInt if result is large, otherwise convert to input type
    if abs(result) > typemax(Int) || T == BigInt
        return BigInt(result)
    else
        return T(result)
    end
end

"""
Get the coefficient of a polynomial at a specific degree.
Returns zero if the degree doesn't exist in the polynomial.
"""
function get_coeff(p::P, deg::D)::C where {C, D, P <: Polynomial{C, D}}
    for t in p
        if t.degree == deg
            return t.coeff
        end
    end
    return zero(C)
end

"""
Chinese Remainder Theorem for TWO polynomials.

Input: [a, b] and [n, m] 
       such that a in Z_n[x], b in Z_m[x], and gcd(n,m)=1.
Output: c in Z_(n*m)[x] 
       such that c == a (mod n) and c == b (mod m).

Algorithm:
    c <- 0
    for k from 0 to max(deg a, deg b) do
        ak <- coefficient of x^k in a (or 0 if k > deg a)
        bk <- coefficient of x^k in b (or 0 if k > deg b)
        ck <- int_crt([ak, bk], [n, m])
        c  <- c + ck * x^k
    return c
"""
function poly_crt_two(a::P, b::P, n::Int, m::Int)::P where {C <: Integer, D, P <: Polynomial{C, D}}
    max_deg = max(degree(a), degree(b))
    
    result_terms = Term{C, D}[]
    
    for k in 0:max_deg
        ak = get_coeff(a, D(k))
        bk = get_coeff(b, D(k))
        
        # Apply integer CRT to coefficients
        ck = int_crt([BigInt(ak), BigInt(bk)], [BigInt(n), BigInt(m)])
        
        if !iszero(ck)
            push!(result_terms, Term(C(ck), D(k)))
        end
    end
    
    return P(result_terms)
end

"""
Chinese Remainder Theorem for multiple polynomials.

Apply CRT iteratively using poly_crt_two.

Input: List of polynomials [p1, p2, ..., pk] and coprime primes [m1, m2, ..., mk]
       where pi is in Z_mi[x]
Output: p in Z_(m1*m2*...*mk)[x] such that p ≡ pi (mod mi) for all i
"""
function poly_crt(polys::Vector{P}, primes::Vector{Int})::P where {C <: Integer, D, P <: Polynomial{C, D}}
    length(polys) == 0 && return P()
    @assert length(polys) == length(primes)
    
    if length(polys) == 1
        return polys[1]
    end
    
    # Start with first polynomial
    result = polys[1]
    M = primes[1]
    
    # Iteratively combine with remaining polynomials
    for i in 2:length(polys)
        result = poly_crt_two(result, polys[i], M, primes[i])
        M *= primes[i]
    end
    
    return result
end

"""
Compute the height (infinity norm) of a polynomial.
height(a) := max{|a0|, |a1|, ..., |an|}
"""
function poly_height(p::P)::BigInt where {P <: Polynomial}
    iszero(p) && return BigInt(0)
    return maximum(abs(BigInt(t.coeff)) for t in p)
end

"""
Symmetric modular reduction for polynomials.
Maps coefficients to range [-m/2, m/2] instead of [0, m-1].

smod(a, m) = (a mod m) if (a mod m) <= m/2 else (a mod m) - m
"""
function smod(p::P, prime::Integer) where {P <: Polynomial}
    p = mod(p, prime)
    s = c -> c > prime ÷ 2 ? c - prime : c
    P(map(t -> Term(s(t.coeff), t.degree), p))
end

"""
Find the next prime number after n.
"""
function next_prime(n::Int)::Int
    candidate = n + 1
    while !is_prime(candidate)
        candidate += 1
    end
    return candidate
end

"""
Simple primality test for small primes.
"""
function is_prime(n::Int)::Bool
    n <= 1 && return false
    n <= 3 && return true
    (n % 2 == 0 || n % 3 == 0) && return false
    
    i = 5
    while i * i <= n
        if n % i == 0 || n % (i + 2) == 0
            return false
        end
        i += 6
    end
    return true
end

"""
Factor a polynomial over Z[x] using CRT-based reconstruction.

This is the main factorization algorithm for Task 7.

Algorithm outline:
1. Choose good primes p1, p2, ... such that f mod pi is square-free with same degree
2. Factor f modulo each prime to get factors in Zp[x]
3. Use CRT to reconstruct factors in Z[x]
4. Verify reconstructed factors divide the original polynomial

Select primes until their product M > 2 * height(f)^deg(f), which provides
sufficient modulus to uniquely reconstruct the integer coefficients.
"""
function factor_crt(f::P)::Vector{Tuple{P, Int}} where {C <: Integer, D, P <: Polynomial{C, D}}
    # Handle simple cases
    degree(f) <= 1 && return [(f, 1)]
    iszero(f) && return [(f, 1)]
    
    # Make f primitive (remove content)
    content_f = gcd(coeffs(f))
    if content_f != 1
        f = P(map(t -> Term(div(t.coeff, content_f), t.degree), f))
    end
    
    # Compute simple bound for coefficient reconstruction
    n = degree(f)
    height_f = poly_height(f)
    
    # Simple bound: M > 2 * height(f)^deg(f)
    # This ensures unique reconstruction of coefficients
    B = 2 * (height_f ^ n)
    
    # Select primes and accumulate until product > B
    primes = Int[]
    M = BigInt(1)
    p = 2
    
    max_attempts = 1000
    attempts = 0
    
    while M <= B && attempts < max_attempts
        p = next_prime(p)
        attempts += 1
        
        # Check if p is a "good" prime (doesn't divide leading coefficient)
        if mod(leading(f).coeff, p) != 0
            # Check if f mod p is square-free and preserves degree
            f_mod_p = mod(f, p)
            if degree(f_mod_p) == degree(f)
                sf = square_free_mod_p(f, p)
                if degree(sf) == degree(f_mod_p)
                    push!(primes, p)
                    M *= p
                end
            end
        end
    end
    
    if isempty(primes)
        # Couldn't find good primes, return original polynomial
        return [(f, 1)]
    end
    
    # Factor modulo each prime
    factorizations = []
    
    for p in primes
        factors_p = factor_mod_p(f, p)
        push!(factorizations, factors_p)
    end
    
    # Use the factorization from the first prime as a guide
    base_factors = factorizations[1]
    
    # Try to reconstruct factors using CRT
    result = Tuple{P, Int}[]
    
    for (base_factor, mult) in base_factors
        # Skip constant factors
        if degree(base_factor) == 0
            continue
        end
        
        # Collect this factor modulo each prime
        factor_images = P[]
        factor_primes = Int[]
        
        for (p, factors_p) in zip(primes, factorizations)
            # Find corresponding factor mod p
            # This is simplified - a full implementation would need sophisticated matching
            for (g_p, m) in factors_p
                if degree(g_p) == degree(base_factor)
                    push!(factor_images, g_p)
                    push!(factor_primes, p)
                    break
                end
            end
        end
        
        if length(factor_images) == length(primes)
            # Reconstruct factor using CRT
            reconstructed = poly_crt(factor_images, factor_primes)
            
            # Apply symmetric mod
            recon_terms = Term{C, D}[]
            half_M = M ÷ 2
            for t in reconstructed
                coeff = t.coeff
                if coeff > half_M
                    coeff = coeff - M
                end
                if !iszero(coeff)
                    push!(recon_terms, Term(coeff, t.degree))
                end
            end
            reconstructed = P(recon_terms)
 
            # Make primitive
            reconstructed = prim_part(reconstructed)
            
            # Verify it divides f
            q, r = pseudo_div_rem(f, reconstructed)
            if iszero(r)
                push!(result, (reconstructed, mult))
            end
        end
    end
    
    # Add back content if needed
    if content_f != 1
        push!(result, (P(Term(C(content_f), zero(D))), 1))
    end
    
    # If no factors found, return original
    if isempty(result)
        result = [(f, 1)]
    end
    
    return result
end

"""
Pseudo-division for polynomials over Z[x] (helper for factor_crt).

Given f, g in Z[x] with g ≠ 0, computes q, r such that:
    lc(g)^(deg(f) - deg(g) + 1) * f = q * g + r
where deg(r) < deg(g) or r = 0.
"""
function pseudo_div_rem(f::P, g::P)::Tuple{P, P} where {C <: Integer, D, P <: Polynomial{C, D}}
    iszero(g) && error("Cannot divide by zero polynomial")
    
    degree(f) < degree(g) && return (zero(P), f)
    
    q = zero(P)
    r = deepcopy(f)
    
    d = degree(f) - degree(g) + 1
    lc_g = leading(g).coeff
    
    # Multiply r by lc(g)^d
    r = r * (lc_g^d)
    
    while degree(r) >= degree(g) && !iszero(r)
        # Compute the leading term of quotient
        lc_r = leading(r).coeff
        deg_r = degree(r)
        deg_g = degree(g)
        
        # Create quotient term
        q_term = Term(lc_r, deg_r - deg_g)
        q = q + P(q_term)
        
        # Subtract from remainder
        r = r - P(q_term) * g
        
        # Trim to remove zero leading terms
        r = trim!(r)
    end
    
    return (q, r)
end
