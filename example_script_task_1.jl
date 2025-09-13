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
