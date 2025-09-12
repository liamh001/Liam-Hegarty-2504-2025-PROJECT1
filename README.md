# 2504_2025_project1 - Polynomial Factorization

*Note: if you have issues with the below commands see the final section on [getting started with the repository](#faq-for-getting-started).*

*Note: `]` means enter `pkg` mode by pressing `]` once. You do **not** need to press it repeatedly.*

This project implements polynomial arithmetic and polynomial factorization for polynomials with integer coefficients.

Students are supposed to create a mirror of the project and create their modifications and improvements according Project1 description. This repository is similar to repositories of previous years, yet has some differences.

To load all functionality, in the directory of the repo:

```
] activate .
```

```
] instantiate
```

```
julia> include("poly_factorization_project.jl")
```

You may then use functionality such as,

```
julia> gcd_mod_p(rand(PolynomialDense) + rand(PolynomialDense), rand(PolynomialDense), 101)
```

To execute all unit tests run:

```
julia> include("test/runtests.jl")
```

You may see examples in `example_script.jl` and run that script line by line.

## FAQ for Getting Started

If you have issues with the above instructions regarding setting up your environment (in particular, `activate`, `instantiate` or `include...`), try the following:

* If `] activate .` worked, try running `] up` (this will ensure the dependencies are versioned correctly for you).

* If `] activate .` did not work, try running `] active raw"PATH"` where `PATH` is the full filepath to the project directory.

* You may also want to try the commands, `] resolve`, `] update` and then again `] instantiate` later.
