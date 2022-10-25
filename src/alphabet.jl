
struct Alphabet
    name :: String
    chars :: Vector{Char}
    wildcard_chars :: Vector{Char}
end

Alphabet(name::AbstractString, chars::AbstractVector{Char}; wildcard_chars="") =
    Alphabet(name, chars, collect(wildcard_chars))
Alphabet(name::AbstractString, chars::AbstractString; wildcard_chars="") =
    Alphabet(name, collect(chars), collect(wildcard_chars))
Alphabet(chars::AbstractString; wildcard_chars="") =
    Alphabet("", collect(chars); wildcard_chars)
Alphabet(chars::AbstractVector{Char}; wildcard_chars="") =
    Alphabet("", chars; wildcard_chars)

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
