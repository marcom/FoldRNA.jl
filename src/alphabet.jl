
struct Alphabet
    name :: String
    chars :: String
end

Alphabet(chars) = Alphabet("", chars)

Base.length(a::Alphabet) = length(a.chars)

encode(a::Alphabet, s) = [findfirst(c, a.chars) for c in s] :: Vector{Int}

decode(a::Alphabet, iter) = join([a.chars[i] for i in iter])
