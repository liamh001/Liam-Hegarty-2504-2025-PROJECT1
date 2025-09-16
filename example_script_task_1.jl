"""
using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

Pkg.instantiate();




######## Task 1.1
x = x_poly(PolynomialDense)

f = x^4 + 4x^2 +(-1)x + 11

g = x^6 + x +(-1)

h = 34x^2 +(-1)x +(-1)

#1.2
sum_fg = +(f,g)
prod_fh = *(f,h)
prod_gh = *(g,h)

#1.3
#product rule manually
product_rule = +(*(f, derivative(g)), *(g,derivative(f)))

@assert ==(derivative(*(f,g)), product_rule)

#1.4
h_inverse(p::Int)  = 
"""

using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

Pkg.instantiate();

# ----------Task 1.1-----------
x = x_poly(PolynomialDense)
f = x^4 + 4x^2 + (-1)x + 11
g = x^6 + x + (-1)
h = 34x^2 + (-1)x + (-1)

# -----------Task 1.2----------
sum_fg = f + g
prod_fh = f * h
prod_gh = g * h

println("f + g = ", sum_fg)
println("f * h = ", prod_fh)
println("g * h = ", prod_gh)
println()


#---------- Task 1.3-------------
fg_product = f * g
derivative_fg = derivative(fg_product)

# Manually
f_prime = derivative(f)
g_prime = derivative(g)
product_rule_result = f_prime * g + f * g_prime

println("Are they equal? ", derivative_fg == product_rule_result)

# --------- Task 1.4
h_inverse(p::Int) = x -> int_inverse_mod(x, p)

function manual_poly_mod(poly::PolynomialDense, n::Integer)
    result_terms = Term{Int,Int}[]
    
    for term in poly
        new_coeff = term.coeff % n
        
        if new_coeff < 0
            new_coeff = new_coeff + n
        end
        
        if new_coeff != 0
            push!(result_terms, Term(new_coeff, term.degree))
        end
    end
    
    # Create new polynomial from modified terms
    return PolynomialDense(result_terms)
end



# Test the implementation
inv_mod_17 = h_inverse(17)

# Compute some inverse values
inv_5 = inv_mod_17(5)   # Should be 7
inv_3 = inv_mod_17(3)   # Should be 6
inv_10 = inv_mod_17(10) # Should be 12

# Test some values
println("Testing h_inverse with prime 17:")
println("Inverse of 5 mod 17: ", inv_5)
println("Inverse of 3 mod 17: ", inv_3)
println("Inverse of 10 mod 17: ", inv_10)

# Verify the inverse property
check_5 = mod(5 * inv_5, 17)
check_3 = mod(3 * inv_3, 17)
check_10 = mod(10 * inv_10, 17)

println("\nVerifying inverse property:")
println("5 * ", inv_5, " mod 17 = ", check_5)  # Should be 1
println("3 * ", inv_3, " mod 17 = ", check_3)  # Should be 1
println("10 * ", inv_10, " mod 17 = ", check_10)  

# --------------------Task 1.5------------------------
println("=== Task 1.5: Modular GCD ===")
fh = f * h
gh = g * h

for p in [5, 11, 13]
    gcd_result = gcd_mod_p(fh, gh, p)
    println("gcd_mod_p(f*h, g*h, ", p, ") = ", gcd_result)
end
println()



# Pretty Printing
function pretty_print(p::PolynomialDense)
        result = ""
           
        for i in length(p.terms):-1:1
            coeff = p.terms[i].coeff
            coeff == 0 && continue
               
            deg = i - 1
               
            # Term string
            term = abs(coeff) == 1 && deg > 0 ? "" : string(abs(coeff))
            term *= deg == 0 ? "" : deg == 1 ? "x" : "x^$deg"
               
            # Concat with result
            if isempty(result)
                result = coeff < 0 ? "-$term" : term
            else
                result *= coeff < 0 ? " - $term" : " + $term"
            end
        end
           
        println(isempty(result) ? "0" : result)
end

        
        
            
        
