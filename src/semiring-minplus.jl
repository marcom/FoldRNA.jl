# (min, +) semiring

export MinPlusSR

struct MinPlusSR{T<:Number} <: Number
    val :: T
end

Base.:(==)(x::MinPlusSR{T}, y::MinPlusSR{T}) where {T} = x.val == y.val
Base.:+(x::MinPlusSR{T}, y::MinPlusSR{T}) where {T} = MinPlusSR{T}(min(x.val, y.val))
Base.:*(x::MinPlusSR{T}, y::MinPlusSR{T}) where {T} = MinPlusSR{T}(x.val + y.val)
Base.zero(::Type{MinPlusSR{T}}) where {T} = MinPlusSR{T}(typemax(T))
# TODO: for integers, we can't use typemax(T) as additions (mul in semiring) would overflow
Base.zero(::Type{MinPlusSR{T}}) where {T<:Integer} = MinPlusSR{T}(T(round(typemax(T)/8)))
Base.one(::Type{MinPlusSR{T}}) where {T} = MinPlusSR{T}(zero(T))
Base.:+(x::MinPlusSR{T}, y::Number) where {T} = x + MinPlusSR{T}(T(y))
Base.:+(x::Number, y::MinPlusSR{T}) where {T} = MinPlusSR{T}(T(x)) + y
Base.:*(x::MinPlusSR{T}, y::Number) where {T} = x * MinPlusSR{T}(T(y))
Base.:*(x::Number, y::MinPlusSR{T}) where {T} = MinPlusSR{T}(T(x)) * y
Base.show(io::IO, ::MIME"text/plain", x::MinPlusSR{T}) where {T} = print(io, "$(x.val)")
