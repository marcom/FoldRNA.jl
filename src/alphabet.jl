
struct Alphabet
    name :: String
    chars :: Vector{Char}
end

Alphabet(name::AbstractString, chars::AbstractString) = Alphabet(name, collect(chars))
Alphabet(chars::AbstractString) = Alphabet("", collect(chars))
Alphabet(chars::AbstractVector{Char}) = Alphabet("", chars)

Base.length(a::Alphabet) = length(a.chars)
Base.:(==)(a::Alphabet, b::Alphabet) = (a.name == b.name && a.chars == b.chars)

encode(a::Alphabet, s) = [findfirst(isequal(c), a.chars) for c in s] :: Vector{Int}
encode(a::Alphabet, c::AbstractChar) = encode(a, string(c))
decode(a::Alphabet, iter) = join(a.chars[i] for i in iter)
