export Pairtable, StrandInfo

const UNPAIRED = 0
const UNPAIRED_CHAR = '.'
const BRACKET_OPEN  = "([{<"
const BRACKET_CLOSE = ")]}>"
const NICK_CHAR = '+'

Base.@kwdef struct StrandInfo
    startidx   :: Int
    endidx     :: Int
    iscircular :: Bool = false
end

struct Pairtable
    pairs   :: Vector{Int}
    strands :: Vector{StrandInfo}
end

function Pairtable(strandlengths::Integer...)
    # TODO: strands default to non-circular
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
        strands[i] = StrandInfo(; startidx, endidx, iscircular=false)
    end
    return Pairtable(pairs, strands)
end

Base.length(pt::Pairtable) = length(pt.pairs)
Base.:(==)(a::Pairtable, b::Pairtable) = (a.pairs == b.pairs) && (a.strands == b.strands)

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
function Pairtable(dbn::String)
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
