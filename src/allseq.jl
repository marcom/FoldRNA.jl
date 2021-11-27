export allseq

struct AllseqIter
    pt :: Pairtable
    bases
    basepairs
    nbases :: Int
    nbasepairs :: Int
end

# TODO: minor inefficiency: collect(allseq(n)) calculates
#       the number of unpaired and basepaired positions twice
AllseqIter(pt, bases, basepairs) =
    AllseqIter(pt, bases, collect(basepairs), length(bases), length(basepairs))

Base.length(it::AllseqIter) = numseq(it.pt; nbases=it.nbases, nbasepairs=it.nbasepairs)

Base.eltype(::AllseqIter) = String

function decode(it::AllseqIter, counts_up, counts_bp)
    n = length(it.pt)
    seq = Vector{Char}(undef, n)
    i_up, i_bp = 1, 1
    for i = 1:n
        if isunpaired(it.pt, i)
            seq[i] = it.bases[counts_up[i_up]]
            i_up += 1
        elseif isbpopening(it.pt, i)
            bp1, bp2 = it.basepairs[counts_bp[i_bp]]
            seq[i] = bp1
            seq[it.pt.pairs[i]] = bp2
            i_bp += 1
        end
    end
    return join(seq)
end

function Base.iterate(it::AllseqIter)
    length(it.pt) == 0 && return nothing
    counts_up = ones(Int, numunpaired(it.pt))
    counts_bp = ones(Int, numbasepaired(it.pt))
    seq = decode(it, counts_up, counts_bp)
    state = counts_up, counts_bp
    return seq, state
end

function Base.iterate(it::AllseqIter, state)
    counts_up, counts_bp = state
    # increment counts
    i, i_up, i_bp = length(it.pt), length(counts_up), length(counts_bp)
    while i >= 1
        if isunpaired(it.pt, i)
            if counts_up[i_up] < it.nbases
                counts_up[i_up] += 1
                break
            end
            counts_up[i_up] = 1
            i_up -= 1
        elseif isbpopening(it.pt, i)
            if counts_bp[i_bp] < it.nbasepairs
                counts_bp[i_bp] += 1
                break
            end
            counts_bp[i_bp] = 1
            i_bp -= 1
        end
        i -= 1
    end
    i < 1 && return nothing
    seq = decode(it, counts_up, counts_bp)
    state = counts_up, counts_bp
    return seq, state
end

"""
    allseq(pt::Pairtable; [bases], [basepairs])
    allseq(dbn::String; [bases], [basepairs])

Iterator over all sequences compatible with a secondary structure
given as a Pairtable `pt` or a string `dbn` in dot-bracket notation.
Optionally, `bases` and `basepairs` to be used can be specified.
"""
allseq(pt::Pairtable; bases=DEFAULT_BASES, basepairs=DEFAULT_BASEPAIRS) =
    AllseqIter(pt, bases, basepairs)
allseq(dbn::String; bases=DEFAULT_BASES, basepairs=DEFAULT_BASEPAIRS) =
    allseq(Pairtable(dbn); bases, basepairs)

"""
    allseq(n::Integer; [bases])

Iterator over all sequences of length `n`.  Optionally, the `bases` to
be used can be specified.
"""
allseq(n::Integer; bases=DEFAULT_BASES) = allseq(Pairtable("."^n); bases)

