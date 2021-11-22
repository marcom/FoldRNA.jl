export Pairtable, StrandInfo

Base.@kwdef struct StrandInfo
    startidx   :: Int
    endidx     :: Int
    iscircular :: Bool = false
end

struct Pairtable
    pairs   :: Vector{Int}
    strands :: Vector{StrandInfo}
end

const UNPAIRED = 0

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
