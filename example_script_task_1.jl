
using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

Pkg.instantiate();

# ----------Task 1.1-----------
x = x_poly(PolynomialDense)
f = x^4 + 4x^2 + (-1)x + 11
g = x^6 + x + (-1)
h = 32x^2 + (-1)x + (-1)

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

# --------- Task 1.4 --------------
product = f * h
for p in [5, 17, 101]
    quotient = div_mod_p(product, h, p)
    
    f_mod_p = mod(f, p)
    
    println("\nModulo $p:")
    println("  (f * h) ÷ h mod $p = ", quotient)
    println("  f mod $p = ", f_mod_p)
    println("  Equal? ", quotient == f_mod_p)
end

# --------------------Task 1.5------------------------
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
               
            # Term str
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



        
        
            
        
