export Fold, energy, mfe, partfn, bpp, prob_of_struct

# Fold

struct Fold{M}
    # TODO: circular strand, multiple strands
    seq :: Vector{Int}
    model :: M
end

Fold(seq::AbstractString, model=RNA_TURNER2004) = Fold(encode(model, seq), model)

Base.length(fold::Fold) = length(fold.seq)

function Base.show(io::IO, mime::MIME"text/plain", fold::Fold)
    println(io, "Fold: $(length(fold)) nts, 1 strand, $(fold.model.name)")
    print(  io, " Seq: ", decode(fold.model, fold.seq))
end

# Note: only implemented for LoopModel at the moment
bptype(fold::Fold{M}, i::Integer, j::Integer) where {M <: LoopModel} =
    bptype(fold.model, fold.seq[i], fold.seq[j])

canbp(fold::Fold{M}, i, j) where {M <: LoopModel} = bptype(fold, i, j) != 0


# functions for ::Fold{BpModel}
# energy, mfe, partfn, bpp, prob_of_struct

function energy(fold::Fold{M}, pt::Pairtable) where {T, M <: BpModel{T}}
    # TODO: should we check here if any base-pair distance is smaller
    #       than hpmin?
    n = length(fold)
    n == length(pt) ||
        throw(ArgumentError("Fold and Pairtable must have same length"))
    en = zero(T)
    for i = 1:n
        !isbpopening(pt, i) && continue
        en += score(fold, Basepair(i, pt.pairs[i]))
    end
    return en * fold.model.unit
end

energy(fold::Fold, dbn::AbstractString) = energy(fold, Pairtable(dbn))

# TODO: backtrack missing
function mfe(fold::Fold{M}) where {T, M <: BpModel{T}}
    A = bpmodel(MinPlusSR{T}, fold; fold.model.hpmin,
                bp = (f, i, j) -> score(f, Basepair(i, j))
        )
    en_mfe = A[1, length(fold)].val
    return fold.model.unit * en_mfe
end

function partfn(fold::Fold{M}) where {T, M <: BpModel{T}}
    A = bpmodel(LogSR{T}, fold; fold.model.hpmin,
                bp = (f, i, j) -> score_exp(f, Basepair(i,j))
        )
    logQ = A[1, length(fold)].val
    return - fold.model.RT * logQ
end

function bpp(fold::Fold{M}) where {T, M <: BpModel{T}}
    A = bpmodel(LogSR{T}, fold; fold.model.hpmin,
                bp = (f, i, j) -> score_exp(f, Basepair(i,j))
        )
    p = bpmodel_bpp(A, fold; hpmin=fold.model.hpmin,
                    bp = (f, i, j) -> score_exp(f, Basepair(i, j)))
    return p
end

function prob_of_struct(fold::Fold, pt::Pairtable)
    RT = fold.model.RT
    en = energy(fold, pt)
    logQ = -partfn(fold) / RT
    logp = -en/RT - logQ
    return exp(logp)
end

prob_of_struct(fold::Fold, dbn::AbstractString) =
    prob_of_struct(fold, Pairtable(dbn))


# functions for ::Fold{LoopModel}
# energy, mfe, partfn, bpp, prob_of_struct

# TODO: backtrack missing
function mfe(fold::Fold{M}) where {T, M <: LoopModel{T}}
    A, Ab, Am, Am1 = loopmodel(MinPlusSR{T}, fold; hpmin=fold.model.hpmin, maxintloop=-1)
    en_mfe = A[1, length(fold)].val
    return fold.model.unit * en_mfe
end

function partfn(fold::Fold{M}) where {T, M <: LoopModel{T}}
    A, Ab, Am, Am1 = loopmodel(LogSR{T}, fold; hpmin=fold.model.hpmin, maxintloop=-1)
    logQ = A[1, length(fold)].val
    return - fold.model.RT * logQ
end
