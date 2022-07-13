
struct Alphabet
    name :: String
    chars :: Vector{Char}
end

Alphabet(name::AbstractString, chars::AbstractString) = Alphabet(name, collect(chars))
Alphabet(chars::AbstractString) = Alphabet("", collect(chars))
Alphabet(chars::AbstractVector{Char}) = Alphabet("", chars)

Base.length(a::Alphabet) = length(a.chars)
Base.:(==)(a::Alphabet, b::Alphabet) = (a.name == b.name && a.chars == b.chars)

function encode(a::Alphabet, s)
    res = [findfirst(isequal(c), a.chars) for c in s]
    length(res) == 0 && return Int[]
    if any(i -> isnothing(i), res)
        throw(ArgumentError("some input character(s) not found in alphabet"))
    end
    return res :: Vector{Int}
end
encode(a::Alphabet, c::AbstractChar) = encode(a, string(c))
decode(a::Alphabet, iter) = join(a.chars[i] for i in iter)
