
# Fold

struct Fold{M}
    # TODO: circular strand, multiple strands
    seq :: Vector{Int}
    model :: M
end

Fold(seq::AbstractString, model) = Fold(encode(model, seq), model)

Base.length(fold::Fold) = length(fold.seq)

function Base.show(io::IO, mime::MIME"text/plain", fold::Fold)
    println(io, "Fold: $(length(fold)) nts, 1 strand, $(fold.model.name)")
    print(  io, " Seq: ", decode(fold.model, fold.seq))
end


# score, score_exp for Fold

struct Basepair
    i :: Int
    j :: Int
end

function score(f::Fold{M}, bp::Basepair) where {M <: BpModel}
    si = f.seq[bp.i]
    sj = f.seq[bp.j]
    BP = f.model.score_bp
    return BP[si, sj]
end

function score_exp(f::Fold{M}, bp::Basepair) where {M <: BpModel}
    si = f.seq[bp.i]
    sj = f.seq[bp.j]
    BP = f.model.score_bp
    RT = ustrip(f.model.RT)
    return exp(- BP[si, sj] / RT)
end


# energy, partfn for Fold

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

function partfn(fold::Fold{M}) where {T, M <: BpModel{T}}
    A = bpmodel(LogSR{T}, fold; fold.model.hpmin,
                bp = (f, i, j) -> score_exp(f, Basepair(i,j))
        )
    logQ = A[1, length(fold)].val
    return - fold.model.RT * logQ
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

