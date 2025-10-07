#############################################################################
# Conversion functions between Integer and ZModP polynomial coefficients
#############################################################################

"""
Convert a polynomial with integer coefficients to a polynomial with ZModP coefficients
"""
function to_zmodp_poly(p::P, prime::Integer) where {C <: Integer, D, P <: Polynomial{C, D}}
    PT = typeof(p)
    terms = Term{ZModP{C, prime}, D}[]
    
    for t in p
        if !iszero(t) && t.degree <= degree(p)  # Filter orphaned terms
            push!(terms, Term(ZModP{C, prime}(t.coeff), t.degree))
        end
    end
    
    # Determine the correct polynomial type
    if PT <: PolynomialDense
        result = PolynomialDense{ZModP{C, prime}, D}(terms)
    elseif PT <: PolynomialSparse
        result = PolynomialSparse{ZModP{C, prime}, D}(terms)
    else
        error("Unknown polynomial type")
    end
    return trim!(result)
end

# Don't convert if already ZModP - just return as-is
function to_zmodp_poly(p::P, prime::Integer) where {T <: Integer, N, D, P <: Polynomial{ZModP{T, N}, D}}
    return p
end

"""
Convert a polynomial with ZModP coefficients back to a polynomial with integer coefficients
"""
function from_zmodp_poly(p::P, ::Type{C}) where {T <: Integer, N, D, P <: Polynomial{ZModP{T, N}, D}, C <: Integer}
    PT = typeof(p)
    terms = Term{C, D}[]
    
    for t in p
        if !iszero(t) && t.degree <= degree(p)  # Filter orphaned terms
            push!(terms, Term(C(t.coeff.val), t.degree))
        end
    end
    
    # Determine the correct polynomial type
    if PT <: PolynomialDense
        result = PolynomialDense{C, D}(terms)
    elseif PT <: PolynomialSparse
        result = PolynomialSparse{C, D}(terms)
    else
        error("Unknown polynomial type")
    end
    return trim!(result)
end

#############################################################################
# Modified functions for ZModP{T, N} coefficient types
#############################################################################

"""
Returns the factors of `f` in an array of tuples. 

For a given tuple (g, n) in the returned array, g is an irreducible factor
of f and the multiplicity of g in f is n.

The returned array may also contain a constant for reconstruction of the
leading coefficient.
"""
function factor(::Type{ZModP{T, N}}, f::P)::Vector{Tuple{P, Int}} where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    # Handle trivial cases
    degree(f) <= 1 && return [(f, 1)]
    
    ret_val = Tuple{P, Int}[]
    
    # IMPORTANT: Store the leading coefficient of the ORIGINAL polynomial
    original_leading_coeff = leading(f).coeff
    
    # Make f square-free
    sqr_fr_poly = square_free(C, f)
    
    # Make the square-free polynomial monic
    fp = sqr_fr_poly
    sqr_free_leading_coeff = leading(fp).coeff
    if sqr_free_leading_coeff != one(C)
        fp = fp * inv(sqr_free_leading_coeff)
    end
    
    # Distinct degree factorization
    dds = dd_factor(C, fp)
    
    # Split each distinct degree component
    for (k, dd) in enumerate(dds)
        if !iszero(dd) && degree(dd) > 0
            sp = dd_split(C, dd, k)
            for mp in sp
                if !iszero(mp) && degree(mp) > 0
                    # Make monic
                    if leading(mp).coeff != one(C)
                        mp = mp * inv(leading(mp).coeff)
                    end
                    push!(ret_val, (mp, multiplicity(C, f, mp)))
                end
            end
        end
    end
    
    # Append the leading coefficient of the ORIGINAL polynomial
    if original_leading_coeff != one(C)
        push!(ret_val, (P(Term(original_leading_coeff, zero(D))), 1))
    end
    
    return ret_val
end

"""
Returns the quotient and remainder of num divided by den (where num/den have the 
same concrete type). 
"""
function div_rem(::Type{ZModP{T, N}}, num::P, den::P)::Tuple{P, P} where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    iszero(den) && throw(DivideError())
    iszero(num) && return zero(P), zero(P)
    
    q = P()
    r = deepcopy(num)
    
    while degree(r) >= degree(den) && !iszero(r)
        h_term = div(leading(r), leading(den))
        h = P(h_term)
        
        r = r - h * den
        q = q + h
        
        trim!(r)
        trim!(q)
    end
    
    return trim!(q), trim!(r)
end

"""
Distinct degree factorization.

Given a square free polynomial `f` returns a list, `g` such that `g[k]` is a product 
of irreducible polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2, such that the 
product of the list is equal to `f`.
"""
function dd_factor(::Type{ZModP{T, N}}, f::P)::Array{P} where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    x = x_poly(P)
    w = deepcopy(x)
    g = Array{P}(undef, degree(f))
    
    # Looping over degrees
    for k in 1:degree(f)
        # w = w^p mod f where p is the characteristic N
        w = rem(C, w^N, f)
        g[k] = gcd(C, w - x, f)
        if !iszero(g[k]) && degree(g[k]) > 0
            f = div(C, f, g[k])
        end
    end
    
    # Edge case for final factor
    if !iszero(f) && degree(f) > 0
        push!(g, f)
    end
    
    # Filter out zero polynomials and return
    return filter(!iszero, g)
end

"""
Distinct degree split.

Returns a list of irreducible polynomials of degree `d` so that the product of
that list is the polynomial `f`.
"""
function dd_split(::Type{ZModP{T, N}}, f::P, d::Integer)::Vector{P} where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    # Base cases
    degree(f) == d && return [f]
    degree(f) == 0 && return P[]
    
    # Generate random polynomial
    w = rand(P, degree = degree(f) - 1, monic = false)
    
    # Compute w^((p^d - 1) / 2) mod f
    n_power = (BigInt(N)^d - 1) ÷ 2
    h = rem(C, w^n_power, f)
    
    # Compute gcd(h - 1, f)
    g = gcd(C, h - one(P), f)
    
    # If g is trivial, try again with a different random polynomial
    if iszero(g) || degree(g) == 0 || degree(g) == degree(f)
        return dd_split(C, f, d)
    end
    
    # Recursive split
    g_bar = div(C, f, g)
    return vcat(dd_split(C, g, d), dd_split(C, g_bar, d))
end

""" 
Returns the quotient of num divided by den (where num/den have the 
same concrete type) 
"""
div(::Type{ZModP{T, N}}, num::P, den::P) where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}} = first(div_rem(C, num, den))

""" 
Returns the remainder of num divided by den (where num/den have the 
same concrete type) 
"""
rem(::Type{ZModP{T, N}}, num::P, den::P) where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}} = last(div_rem(C, num, den))

"""
The extended euclid algorithm for polynomials (of the same concrete subtype).
"""
function extended_euclid_alg(::Type{ZModP{T, N}}, f::P, g::P) where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    return ext_euclid_alg(f, g, (a,b) -> rem(C, a, b), (a,b) -> div(C, a, b))
end

"""
The greatest common divisor of two polynomials (of the same concrete subtype).
"""
gcd(::Type{ZModP{T, N}}, f::P, g::P) where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}} = extended_euclid_alg(C, f, g) |> first

"""
Square free factorization for polynomials in Zp[x] (Yun's algorithm)
"""
function square_free(::Type{ZModP{T, N}}, f::P) where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    # Remove minimum degree (in case there's a factor of x)
    min_deg = last(f).degree
    vt = filter(t -> !iszero(t), collect(f))
    f = P(map(t -> Term(t.coeff, t.degree - min_deg), vt))
    
    # Compute the gcd of f and f'
    der_f = derivative(f)
    sqr_part = gcd(C, f, der_f)
    
    # If gcd is trivial, f is already square-free
    iszero(sqr_part) && return f * (min_deg > zero(D) ? x_poly(P) : one(P))
    
    # Remove factors with multiplicity > 1
    sqr_free = div(C, f, sqr_part)
    
    # Add the factor of x back if necessary
    if min_deg > zero(D)
        sqr_free *= x_poly(P)
    end
    
    return sqr_free
end

"""
Compute the number of times g divides f.
"""
function multiplicity(::Type{ZModP{T, N}}, f::P, g::P)::Integer where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    degree(gcd(C, f, g)) == 0 && return 0
    return 1 + multiplicity(C, div(C, f, g), g)
end

#############################################################################
# Refactored _mod_p functions to use ZModP
#############################################################################

"""
Division with remainder modulo a prime - refactored to use ZModP
"""
function div_rem_mod_p(num::P, den::P, prime::Integer)::Tuple{P, P} where {C <: Integer, D, P <: Polynomial{C, D}}
    # Convert to ZModP polynomials
    num_zp = to_zmodp_poly(num, prime)
    den_zp = to_zmodp_poly(den, prime)
    
    # Perform division in Zp[x]
    q_zp, r_zp = div_rem(ZModP{C, prime}, num_zp, den_zp)
    
    # Convert back
    q = from_zmodp_poly(q_zp, C)
    r = from_zmodp_poly(r_zp, C)
    
    return q, r
end

"""
Division with remainder for polynomials that already have ZModP coefficients
"""
function div_rem_mod_p(num::P, den::P, prime::Integer) where {T <: Integer, N, D, P <: Polynomial{ZModP{T, N}, D}}
    # Already in Zp[x], use the ZModP division directly
    return div_rem(ZModP{T, N}, num, den)
end

"""
Extended Euclidean algorithm modulo a prime - refactored to use ZModP
"""
function extended_euclid_alg_mod_p(a::P, b::P, prime::Integer) where {C <: Integer, D, P <: Polynomial{C, D}}
    # Convert to ZModP polynomials
    a_zp = to_zmodp_poly(a, prime)
    b_zp = to_zmodp_poly(b, prime)
    
    # Perform extended Euclidean algorithm in Zp[x]
    g_zp, s_zp, t_zp = extended_euclid_alg(ZModP{C, prime}, a_zp, b_zp)
    
    # Convert back
    g = from_zmodp_poly(g_zp, C)
    s = from_zmodp_poly(s_zp, C)
    t = from_zmodp_poly(t_zp, C)
    
    return g, s, t
end

# If already ZModP, use the ZModP version directly
function extended_euclid_alg_mod_p(a::P, b::P, prime::Integer) where {T <: Integer, N, D, P <: Polynomial{ZModP{T, N}, D}}
    return extended_euclid_alg(ZModP{T, N}, a, b)
end

"""
Square-free factorization modulo a prime - refactored to use ZModP
"""
function square_free_mod_p(f::P, prime::Integer) where {C <: Integer, D, P <: Polynomial{C, D}}
    # Convert to ZModP polynomial
    f_zp = to_zmodp_poly(f, prime)
    
    # Compute square-free in Zp[x]
    sqf_zp = square_free(ZModP{C, prime}, f_zp)
    
    # Convert back
    return from_zmodp_poly(sqf_zp, C)
end

"""
Factorization modulo a prime - refactored to use ZModP
"""
function factor_mod_p(f::P, prime::Integer)::Vector{Tuple{P, Int}} where {C <: Integer, D, P <: Polynomial{C, D}}
    # Convert to ZModP polynomial
    f_zp = to_zmodp_poly(f, prime)
    
    # Factor in Zp[x]
    factors_zp = factor(ZModP{C, prime}, f_zp)
    
    # Convert back
    factors = Tuple{P, Int}[]
    for (g_zp, mult) in factors_zp
        g = from_zmodp_poly(g_zp, C)
        push!(factors, (g, mult))
    end
    
    return factors
end

"""
Multiplicity modulo a prime - refactored to use ZModP
"""
function multiplicity_mod_p(f::P, g::P, prime::Integer)::Int where {C <: Integer, D, P <: Polynomial{C, D}}
    # Convert to ZModP polynomials
    f_zp = to_zmodp_poly(f, prime)
    g_zp = to_zmodp_poly(g, prime)
    
    # Compute multiplicity in Zp[x]
    return multiplicity(ZModP{C, prime}, f_zp, g_zp)
end

"""
Distinct degree factorization modulo a prime - refactored to use ZModP
"""
function dd_factor_mod_p(f::P, prime::Integer)::Array{P} where {C <: Integer, D, P <: Polynomial{C, D}}
    # Convert to ZModP polynomial
    f_zp = to_zmodp_poly(f, prime)
    
    # Perform distinct degree factorization in Zp[x]
    dds_zp = dd_factor(ZModP{C, prime}, f_zp)
    
    # Convert back
    dds = P[]
    for dd_zp in dds_zp
        push!(dds, from_zmodp_poly(dd_zp, C))
    end
    
    return dds
end

"""
Distinct degree split modulo a prime - refactored to use ZModP
"""
function dd_split_mod_p(f::P, d::Integer, prime::Integer)::Vector{P} where {C <: Integer, D, P <: Polynomial{C, D}}
    # Convert to ZModP polynomial
    f_zp = to_zmodp_poly(f, prime)
    
    # Perform distinct degree split in Zp[x]
    splits_zp = dd_split(ZModP{C, prime}, f_zp, d)
    
    # Convert back
    splits = P[]
    for sp_zp in splits_zp
        push!(splits, from_zmodp_poly(sp_zp, C))
    end
    
    return splits
end
