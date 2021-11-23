export Pairtable, StrandInfo, numseq, randseq

const UNPAIRED = 0
const UNPAIRED_CHAR = '.'
const BRACKET_OPEN  = "([{<"
const BRACKET_CLOSE = ")]}>"
const NICK_CHAR = '+'

"""
    StrandInfo(startidx, endidx, iscircular)
    StrandInfo(; startidx, endidx, iscircular)

Strand information used in a Pairtable: start index `startidx` of the
strand in the Pairtable, end index `endidx`, and if the strand
`iscircular`.
"""
Base.@kwdef struct StrandInfo
    startidx   :: Int
    endidx     :: Int
    # TODO: At the moment, strands can be either linear or circular.
    # Should also support lariats.  The end should also be able to
    # link back to any position (not just the first base as in the
    # circular case).
    # Possible solution:
    # linkage_end :: Int   # 0 would mean linear
    iscircular :: Bool = false
end
Base.hash(s::StrandInfo, h::UInt) = hash(s.iscircular, hash(s.endidx, hash(s.startidx, h)))

"""
    StrandInfo(startidx, endidx) = StrandInfo(startidx, endidx, false)

Strand information for a linear (non-circular) strand.
"""
StrandInfo(startidx, endidx) = StrandInfo(startidx, endidx, false)

struct Pairtable
    pairs   :: Vector{Int}
    strands :: Vector{StrandInfo}
end

function Pairtable(strandlengths::Integer...)
    # TODO: Strands default to non-circular.
    nstrand = length(strandlengths)
    if ! (nstrand > 0)
        throw(ArgumentError("No strandlengths given"))
    end
    if ! all(i > 0 for i in strandlengths)
        throw(ArgumentError("All strandlengths must be positive"))
    end
    n = sum(strandlengths)
    pairs = fill(UNPAIRED, n)
    strands = Vector{StrandInfo}(undef, nstrand)
    startidx = 0
    endidx = 0
    for (i, len) in enumerate(strandlengths)
        startidx = endidx + 1
        endidx = endidx + len
        strands[i] = StrandInfo(startidx, endidx)
    end
    return Pairtable(pairs, strands)
end

Base.length(pt::Pairtable) = length(pt.pairs)
Base.:(==)(a::Pairtable, b::Pairtable) = (a.pairs == b.pairs) && (a.strands == b.strands)
Base.hash(a::Pairtable, h::UInt) = hash(a.strands, hash(a.pairs, h))
hasbp(pt::Pairtable, i::Int, j::Int) = pt.pairs[i] == j
isunpaired(pt::Pairtable, i) = pt.pairs[i] == UNPAIRED
isbpopening(pt::Pairtable, i) = pt.pairs[i] != UNPAIRED && i < pt.pairs[i]
isbpclosing(pt::Pairtable, i) = pt.pairs[i] != UNPAIRED && pt.pairs[i] < i
numunpaired(pt::Pairtable) = sum(isunpaired(pt, i) for i = 1:length(pt))
numbasepaired(pt::Pairtable) = sum(isbpopening(pt, i) for i = 1:length(pt))

function Base.String(pt::Pairtable)
    # TODO: pseudoknots not handled or or even detected
    # TODO: circular strands not handled
    n = length(pt)
    s = fill(UNPAIRED_CHAR, n)
    for i in 1:length(pt)
        if isbpopening(pt, i)
            s[i] = BRACKET_OPEN[1]
            s[pt.pairs[i]] = BRACKET_CLOSE[1]
        end
    end
    return join([ String(s[si.startidx:si.endidx]) for si in pt.strands ], NICK_CHAR)
end


"""
    numseq(pt::Pairtable; nbases=DEFAULT_NBASES, nbasepairs=DEFAULT_NBASEPAIRS)
    numseq(dbn::String; nbases=DEFAULT_NBASES, nbasepairs=DEFAULT_NBASEPAIRS)

Number of sequences for a given secondary structure given as a
`Pairtable` or as a string in dot-bracket notation.  Optionally the
number of possible bases `nbases` and basepairs `nbasepairs` can be
specified.
"""
numseq(pt::Pairtable; nbases=DEFAULT_NBASES, nbasepairs=DEFAULT_NBASEPAIRS) =
    big(nbases)^numunpaired(pt) * big(nbasepairs)^numbasepaired(pt)
numseq(dbn::String; nbases=DEFAULT_NBASES, nbasepairs=DEFAULT_NBASEPAIRS) =
    numseq(Pairtable(dbn); nbases, nbasepairs)
numseq(n::Integer; nbases=DEFAULT_NBASES) = big(nbases)^n


"""
    randseq(pt::Pairtable; [bases], [basepairs])
    randseq(dbn::String; [bases], [basepairs])

Random sequence for a secondary structure given as a `Pairtable` or
string in dot-bracket notation.  Optionally, the `bases` and
`basepairs` to be used can be specified.
"""
function randseq(pt::Pairtable;
                 bases=DEFAULT_BASES,
                 basepairs=DEFAULT_BASEPAIRS)
    if length(pt.strands) != 1
        throw(ArgumentError("randseq doesn't handle multiple strands yet"))
    end
    n = length(pt)
    seq = Vector{Char}(undef, n)
    for i = 1:n
        if isunpaired(pt, i)
            seq[i] = rand(bases)
        elseif isbpopening(pt, i)
            a, b = rand(basepairs)
            seq[i] = a
            seq[pt.pairs[i]] = b
        end
    end
    return join(seq)
end
randseq(dbn::AbstractString; bases=DEFAULT_BASES, basepairs=DEFAULT_BASEPAIRS) =
    randseq(Pairtable(dbn); bases, basepairs)

"""
    randseq(n::Integer; [bases])

Random sequence of length `n`.  Optionally, the possible `bases` to be
used can be specified.
"""
randseq(n::Integer; bases=DEFAULT_BASES) = randseq(Pairtable("."^n); bases)


"""
    Pairtable(dbn::String)

Compute the pair table (base pairs and strand extents) from `dbn`, a
nucleic acid secondary structure given in dot-bracket notation.

In dot-bracket notation, unpaired bases are denoted by '.', base pairs
are formed by pairing brackets '()[]{}<>', and strand ends (in a
multi-strand complex) are indicated by '+'.

Throws an exception if the dot-bracket string is malformed (illegal
characters or wrong number of opening or closing parentheses).

# Examples
```jldoctest
julia> Pairtable("(..)")
Pairtable([4, 0, 0, 1], StrandInfo[StrandInfo(1, 4, false)])

julia> Pairtable("(.(.[+.).].)")
Pairtable([11, 0, 7, 0, 9, 0, 3, 0, 5, 0, 1], StrandInfo[StrandInfo(1, 5, false), StrandInfo(6, 11, false)])
```
"""
function Pairtable(dbn::AbstractString)
    # TODO: doesn't support circular strands
    nbracket = length(BRACKET_OPEN)
    nstrand = count(c -> c == NICK_CHAR, dbn) + 1
    n = length(dbn) - (nstrand - 1)
    stack = Vector{Vector{Int}}(undef, nbracket)
    for k = 1:nbracket
        stack[k] = Int[]
    end
    pairs = fill(UNPAIRED, n)
    strands = Vector{StrandInfo}(undef, nstrand)
    for k = 1:nstrand
        strands[k] = StrandInfo(0, 0, false)
    end
    strands[1] = StrandInfo(1, 0, false)

    i = 1
    i_strand = 1
    for c in dbn
        if c == NICK_CHAR
            strands[i_strand] = StrandInfo(strands[i_strand].startidx, i-1, false)
            if strands[i_strand].startidx > strands[i_strand].endidx
                error("complex has a strand with length 0 ending at position $(i + i_strand - 1)")
            end
            i_strand += 1
            strands[i_strand] = StrandInfo(i, 0, false)
            continue
        end
        idx_open = something(findfirst(isequal(c), BRACKET_OPEN), 0)
        idx_close = something(findfirst(isequal(c), BRACKET_CLOSE), 0)
        if c == UNPAIRED_CHAR
            pairs[i] = UNPAIRED
        elseif idx_open > 0
            # opening bracket
            push!(stack[idx_open], i)
        elseif idx_close > 0
            # closing bracket
            if length(stack[idx_close]) > 0
                j = pop!(stack[idx_close])
                pairs[i] = j
                pairs[j] = i
            else
                error("too many closing brackets of type '$(BRACKET_CLOSE[idx_close])' at position $(i + i_strand - 1)")
            end
        else
            error("illegal character '$c' in dot-bracket string")
        end
        i += 1
    end
    for k = 1:nbracket
        if length(stack[k]) != 0
            error("more opening brackets than closing brackets of type '$(BRACKET_OPEN[k])' ($(length(stack[k])) more)")
        end
    end
    strands[nstrand] = StrandInfo(strands[nstrand].startidx, n, false)
    if strands[nstrand].startidx > strands[nstrand].endidx
        error("complex has a strand with length 0 ending at position $n")
    end
    return Pairtable(pairs, strands)
end
