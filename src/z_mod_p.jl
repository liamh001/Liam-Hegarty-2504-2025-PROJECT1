using Primes
import Base: +, -, *, ÷, ^, ==, show, zero, one, inv, abs

struct ZModP{T <: Integer, N} <: Integer
    val::T

    function ZModP{T, N}(a::T) where {T <: Integer, N}
        isprime(N) || error("N = $N must be prime")
        new(mod(a, N))
    end
end

# Outer constructors
ZModP(a::T, ::Val{N}) where {T <: Integer, N} = ZModP{T, N}(a)

# Allow conversion from different integer types
function ZModP{T, N}(a::S) where {T <: Integer, S <: Integer, N}
    ZModP{T, N}(T(a))  # Convert to type T first
end

# Zero and one
zero(::Type{ZModP{T, N}}) where {T <: Integer, N} = ZModP{T, N}(zero(T))
one(::Type{ZModP{T, N}}) where {T <: Integer, N} = ZModP{T, N}(one(T))
zero(a::ZModP{T, N}) where {T <: Integer, N} = zero(ZModP{T, N})
one(a::ZModP{T, N}) where {T <: Integer, N} = one(ZModP{T, N})

# Display
show(io::IO, a::ZModP{T, N}) where {T <: Integer, N} = print(io, a.val)

# Equality
==(a::ZModP{T, N}, b::ZModP{T, N}) where {T <: Integer, N} = a.val == b.val
==(a::ZModP{T, N}, b::S) where {T <: Integer, S <: Integer, N} = a.val == mod(b, N)
==(a::S, b::ZModP{T, N}) where {T <: Integer, S <: Integer, N} = mod(a, N) == b.val

# Addition
+(a::ZModP{T, N}, b::ZModP{T, N}) where {T <: Integer, N} = ZModP{T, N}(a.val + b.val)
+(a::ZModP{T, N}, b::S) where {T <: Integer, S <: Integer, N} = ZModP{T, N}(a.val + b)
+(a::Int, b::ZModP{T, N}) where {T <: Integer, N} = ZModP{T, N}(a + b.val)

# Negation and subtraction
-(a::ZModP{T, N}) where {T <: Integer, N} = ZModP{T, N}(-a.val)
-(a::ZModP{T, N}, b::ZModP{T, N}) where {T <: Integer, N} = ZModP{T, N}(a.val - b.val)
-(a::ZModP{T, N}, b::S) where {T <: Integer, S <: Integer, N} = ZModP{T, N}(a.val - b)
-(a::S, b::ZModP{T, N}) where {T <: Integer, S <: Integer, N} = ZModP{T, N}(a - b.val)

# Multiplication
*(a::ZModP{T, N}, b::ZModP{T, N}) where {T <: Integer, N} = ZModP{T, N}(a.val * b.val)
*(a::ZModP{T, N}, b::S) where {T <: Integer, S <: Integer, N} = ZModP{T, N}(a.val * b)
*(a::S, b::ZModP{T, N}) where {T <: Integer, S <: Integer, N} = ZModP{T, N}(a * b.val)

# Inverse and division
function inv(a::ZModP{T, N})::ZModP{T, N} where {T <: Integer, N}
    a.val == 0 && error("Cannot invert 0 in Z_$N")
    ZModP{T, N}(int_inverse_mod(a.val, N))
end

÷(a::ZModP{T, N}, b::ZModP{T, N}) where {T <: Integer, N} = a * inv(b)
÷(a::ZModP{T, N}, b::S) where {T <: Integer, S <: Integer, N} = a ÷ ZModP{T, N}(b)
÷(a::S, b::ZModP{T, N}) where {T <: Integer, S <: Integer, N} = ZModP{T, N}(a) ÷ b

# Power
function ^(a::ZModP{T, N}, n::S) where {T <: Integer, S <: Integer, N}
    n < 0 && return inv(a)^(-n)
    result = one(ZModP{T, N})
    base = a
    while n > 0
        if n % 2 == 1
            result *= base
        end
        base *= base
        n ÷= 2
    end
    result
end

# Absolute value
abs(a::ZModP{T, N}) where {T <: Integer, N} = a

# Conversion functions
Base.Int(a::ZModP{T, N}) where {T <: Integer, N} = Int(a.val)
Base.BigInt(a::ZModP{T, N}) where {T <: Integer, N} = BigInt(a.val)
