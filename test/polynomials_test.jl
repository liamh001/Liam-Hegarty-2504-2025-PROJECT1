#############################################################################
#############################################################################
#
# This file contains units tests for polynomial operations
#                                                                               
#############################################################################
#############################################################################


"""
Test product of polynomials.
"""
function prod_test_poly(::Type{P};
    N::Int = 3 * 10^2, N_prods::Int = 20, seed::Int = 0
    ) where {C,D, P <: Polynomial{C,D}}
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(P)
        p2 = rand(P)
        prod = p1*p2
        @assert leading(prod) == leading(p1)*leading(p2)
    end
    for _ in 1:N
        p_base = P(Term(one(C),zero(D)))
        for _ in 1:N_prods
            p = rand(P)
            prod = p_base*p
            @assert leading(prod) == leading(p_base)*leading(p)
            p_base = prod
        end
    end
    println("prod_test_poly for $(P) - PASSED")
end

"""
Test derivative of polynomials (as well as product).
"""
function prod_derivative_test_poly(::Type{P};
    N::Int = 10^2,  seed::Int = 0
    ) where {C,D, P <: Polynomial{C,D}}
    Random.seed!(seed)
    for i in 1:N
        p1 = rand(P)
        p2 = rand(P)
        p1d = derivative(p1)
        p2d = derivative(p2)
        
        lhs = p1d * p2 + p1 * p2d
        rhs = derivative(p1 * p2)
        
        if lhs != rhs
            println("FAILED on iteration $i")
            println("p1 = $p1")
            println("p2 = $p2")
            println("lhs = $lhs")
            println("rhs = $rhs")
            
            # ADD THIS DEBUG CODE:
            println("\nLHS terms structure:")
            for (j, t) in enumerate(lhs.terms)
                println("  Position $j (degree $(j-1)): coeff=$(t.coeff), stored_degree=$(t.degree)")
            end
            
            println("\nRHS terms structure:")
            for (j, t) in enumerate(rhs.terms)
                println("  Position $j (degree $(j-1)): coeff=$(t.coeff), stored_degree=$(t.degree)")
            end
            
            println("\nlhs - rhs = $(lhs - rhs)")
            @assert false
        end
    end
    println("prod_derivative_test_poly for $(P) - PASSED")
end
"""
Test division of polynomials modulo p.
"""
function division_test_poly(::Type{P};
    prime::Int = 101, N::Int = 10^4, seed::Int = 0
    ) where {C,D, P <: Polynomial{C,D}}
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(P)
        p2 = rand(P)
        p_prod = p1*p2
        q, r = P(), P()
        try
            q, r = div_rem_mod_p(p_prod, p2, prime)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero( mod(q*p2+r - p_prod, prime) )
    end
    println("division_test_poly for $(P) - PASSED")
end

"""
Test the extended euclid algorithm for polynomials modulo p.
"""
function ext_euclid_test_poly(::Type{P};
    prime::Int=101, N::Int = 10^3, seed::Int = 0
    ) where {C,D, P <: Polynomial{C,D}}
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(P)
        p2 = rand(P)
        g, s, t = extended_euclid_alg_mod_p(p1, p2, prime)
        @assert mod(s*p1 + t*p2 - g, prime) == 0
    end
    println("ext_euclid_test_poly for $(P) - PASSED")
end
