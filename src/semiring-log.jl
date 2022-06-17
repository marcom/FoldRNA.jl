# Log semiring
#
# TODO
# - keep LogSR(x) constructor?
# - `show` LogSR{Float64} as LogSR

export LogSR
using LogExpFunctions: logaddexp

struct LogSR{T<:Number} <: Number
    val :: T
end

#LogSR(x) = LogSR{Float64}(x)
Base.float(x::LogSR{T}) where {T} = exp(x.val)
Base.:(==)(x::LogSR{T}, y::LogSR{T}) where {T} = x.val == y.val
Base.:+(x::LogSR{T}, y::LogSR{T}) where {T} = LogSR{T}(logaddexp(x.val, y.val))
Base.:*(x::LogSR{T}, y::LogSR{T}) where {T} = LogSR{T}(x.val + y.val)
Base.:/(x::LogSR{T}, y::LogSR{T}) where {T} = LogSR{T}(x.val - y.val)
Base.zero(::Type{LogSR{T}}) where {T} = LogSR{T}(typemin(T))
Base.one(::Type{LogSR{T}}) where {T} = LogSR{T}(zero(T))
Base.:+(x::LogSR{T}, y::Number) where {T} = x + LogSR{T}(log(T(y)))
Base.:+(x::Number, y::LogSR{T}) where {T} = LogSR{T}(log(T(x))) + y
Base.:*(x::LogSR{T}, y::Number) where {T} = x * LogSR{T}(log(T(y)))
Base.:*(x::Number, y::LogSR{T}) where {T} = LogSR{T}(log(T(x))) * y
Base.show(io::IO, ::MIME"text/plain", x::LogSR{T}) where {T} = print(io, "exp($(x.val))")
