using LinearAlgebra: dot
using UnicodePlots: heatmap

export FoldPseq

struct FoldPseq{T,M}
    # TODO: circular strand, multiple strands
    pseq :: Matrix{T}
    model :: M
end

FoldPseq(seq::AbstractString, model) = FoldPseq(onehot(model.alphabet, seq), model)

Base.length(fold::FoldPseq) = size(fold.pseq, 2)

function Base.show(io::IO, mime::MIME"text/plain", fold::FoldPseq)
    println(io, "FoldPseq: $(length(fold)) nts, 1 strand, $(fold.model.name)")
    println(io, "      ", String(fold))
    show(io, mime, heatmap(fold.pseq; colormap=:inferno))
    println(io)
    show(io, mime, fold.pseq)
end

String(fold::FoldPseq) = decode(fold.model.alphabet, hotidx(fold.pseq))

function energy(fold::FoldPseq{T,M}, pt::Pairtable) where {T, Tp, M <: BpModel{Tp}}
    # Note: en has to be of type T (coming from FoldPseq) for
    #       forward-mode autodiff to work
    n = length(fold)
    n == length(pt) ||
        throw(ArgumentError("FoldPseq and Pairtable must have same length"))
    en = zero(T)
    for i = 1:n
        !isbpopening(pt, i) && continue
        en += score(fold, Basepair(i, pt.pairs[i]))
    end
    return en * fold.model.unit
end

energy(fold::FoldPseq, dbn::AbstractString) = energy(fold, Pairtable(dbn))

function partfn(fold::FoldPseq{T,M}) where {T, Tp, M <: BpModel{Tp}}
    # Note: LogSR{T} has to use type T (coming from FoldPseq) for
    #       forward-mode autodiff to work
    A = bpmodel(LogSR{T}, fold; fold.model.hpmin,
                bp = (f, i, j) -> score_exp(f, Basepair(i,j))
        )
    logQ = A[1, length(fold)].val
    return - fold.model.RT * logQ
end

function prob_of_struct(fold::FoldPseq, pt::Pairtable)
    RT = fold.model.RT
    en = energy(fold, pt)
    logQ = -partfn(fold) / RT
    logp = -en/RT - logQ
    return exp(logp)
end

prob_of_struct(fold::FoldPseq, dbn::AbstractString) =
    prob_of_struct(fold, Pairtable(dbn))


# onehot, isohot, hotidx for pseq

function onehot(::Type{T}, items, arr) where {T}
    n = length(arr)
    nitem = length(items)
    t = zeros(T, nitem, n)
    for i = 1:n
        k = findfirst(arr[i], items)
        isnothing(k) && error("illegal item $(arr[i]) at position $i")
        t[k, i] = one(T)
    end
    return t
end
onehot(items, arr) = onehot(Float64, items, arr)
onehot(::Type{T}, a::Alphabet, s) where {T} = onehot(T, a.chars, s)
onehot(a::Alphabet, s) = onehot(a.chars, s)

isohot(::Type{T}, d::Integer, n::Integer) where {T} = ones(T, d, n) / d
isohot(d::Integer, n::Integer) = isohot(Float64, d, n)
isohot(::Type{T}, a::Alphabet, n::Integer) where {T} = isohot(T, length(a.chars), n)
isohot(a::Alphabet, n::Integer) = isohot(length(a.chars), n)

# inverse of onehot
# possible alternative names: idxmax (pandas), hotone
function hotidx(pseq::Matrix; dims=1)
    _, cidxs = findmax(pseq; dims)
    otherdim = dims == 1 ? 2 : 1
    return reshape(map(c -> c[dims], cidxs), size(pseq, otherdim))
end
