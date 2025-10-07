#############################################################################
#############################################################################
#
# This file implements factorization 
#                                                                               
#############################################################################
#############################################################################

"""
Factors a polynomial over the field Z_p.

Returns a vector of tuples of (irreducible polynomials (mod p), multiplicity) such that their product of the list (mod p) is f. 
Irreducibles are fixed points on the function `factor`.
"""
#=
function factor_mod_p(f::P, prime::Integer)::Vector{Tuple{P,Integer}} where {P <: Polynomial}
    # Cantor Zassenhaus factorization

    f_modp = mod(f, prime)
    degree(f_modp) ≤ 1 && return [(f_modp,1)]

    ret_val = Tuple{P, Int}[]

    # make f square-free
    sqr_fr_poly = square_free_mod_p(f, prime)

    # Drop to Zp[x]
    fp = mod(sqr_fr_poly, prime)
    # @show "see square free part" fp 

    # make f primitive
    content = mod(gcd(map(t -> t.coeff, fp)), prime)
    fp = mod(fp * int_inverse_mod(content, prime), prime)
    # @show "after square free:", ff

    # make f monic
    old_coeff = leading(fp).coeff
    fp = mod(fp * int_inverse_mod(old_coeff, prime), prime)
    # @show "after monic:", ff

    dds = dd_factor_mod_p(fp, prime)

    for (k,dd) in enumerate(dds)
        sp = dd_split_mod_p(dd, k, prime)
        sp = map((p)->div_mod_p(p, leading(p).coeff, prime), sp) #makes the polynomials inside the list sp, monic
        for mp in sp
            push!(ret_val, (mp, multiplicity_mod_p(f_modp,mp,prime)) )
        end
    end

    #Append the leading coefficient as well
    push!(ret_val, (leading(f_modp).coeff* one(f), 1) )

    return ret_val
end
=#

"""
Expand a factorization.
"""
function expand_factorization(factorization::Vector{Tuple{P, T}})::P where {P <: Polynomial, T <: Integer}
    length(factorization) == 1 && return first(factorization[1])^last(factorization[1])
    return *([first(tt)^last(tt) for tt in factorization]...)
end

"""
Compute the number of times g divides f modulo a prime.
"""
function multiplicity_mod_p(f::P, g::P, prime::Integer)::Integer where {P <: Polynomial}
    degree(gcd_mod_p(f, g, prime)) == 0 && return 0
    return 1 + multiplicity_mod_p(div_mod_p(f, g, prime), g, prime)
end


"""
Distinct degree factorization.

Given a square free polynomial `f` returns a list, `g` such that `g[k]` is a product of irreducible 
polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2, such that the product of the list (mod `prime`) 
is equal to `f` (mod `prime`).
"""
function dd_factor_mod_p(f::P, prime::Integer)::Array{P} where {P <: Polynomial}
    x = x_poly(P)
    w = deepcopy(x)
    g = Array{P}(undef,degree(f)) #Array of polynomials indexed by degree

    #Looping over degrees
    for k in 1:degree(f)
        w = rem_mod_p(pow_mod(w,prime,prime), f, prime)
        g[k] = gcd_mod_p(w - x, f, prime) 
        f = div_mod_p(f, g[k], prime)
    end

    #edge case for final factor
    f != one(P) && push!(g,f)
    
    return g
end

"""
Distinct degree split.

Returns a list of irreducible polynomials of degree `d` so that the product of that list (mod prime) is the polynomial `f`.
"""
function dd_split_mod_p(f::P, d::Integer, prime::Integer)::Vector{P} where {P <: Polynomial}
    f = mod(f,prime)
    degree(f) == d && return [f]
    degree(f) == 0 && return []
    w = rand(P, degree = d, monic = true)
    w = mod(w,prime)
    n_power = (prime^d-1) ÷ 2 # ÷ is integer quotient
    g = gcd_mod_p(pow_mod(w,n_power,prime) - one(f), f, prime)
    g = mod(g, prime)
    # @show g
    ḡ = div_mod_p(f, g, prime) # g\bar + [TAB]
    return vcat(dd_split_mod_p(g, d, prime), dd_split_mod_p(ḡ, d, prime) )
end
