
const DEFAULT_BASES = [
    'A', 'C', 'G', 'U'
]

const DEFAULT_BASEPAIRS = Set([
    ('A', 'U'),
    ('U', 'A'),
    ('C', 'G'),
    ('G', 'C'),
    ('G', 'U'),
    ('U', 'G'),
])

const DEFAULT_NBASES = length(DEFAULT_BASES)
const DEFAULT_NBASEPAIRS = length(DEFAULT_BASEPAIRS)

function default_canbp(seq::AbstractString, i, j)
    a = seq[i]
    b = seq[j]
    return a == 'N' || b == 'N' || (a,b) in DEFAULT_BASEPAIRS
end

function default_canbp_noN(seq::AbstractString, i, j)
    return (seq[i], seq[j]) in DEFAULT_BASEPAIRS
end
