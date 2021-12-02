function exhaustive_mfe(fold::Fold)
    n = length(fold)
    pt_min = Pairtable(n)
    en_min = energy(fold, pt_min)
    seq = decode(fold.model.al, fold.seq)
    for pt in allstruct(seq; fold.model.hpmin)
        en = energy(fold, pt)
        if en < en_min
            en_min = en
            pt_min = deepcopy(pt)
        end
    end
    return en_min, pt_min
end

function exhaustive_partfn(::Type{T}, fold::Fold) where {T}
    RT, hpmin = fold.model.RT, fold.model.hpmin
    seq = decode(fold.model.al, fold.seq)
    Q = zero(T)
    for pt in allstruct(seq; hpmin)
        Q += exp(- energy(fold, pt) / RT)
    end
    return Q
end

function exhaustive_partfn(fold::Fold)
    s = exhaustive_partfn(LogSR{Float64}, fold)
    logQ = s.val
    return - fold.model.RT * logQ
end

function exhaustive_bpp_partfn(::Type{T}, fold::Fold) where {T}
    n = length(fold)
    RT, hpmin = fold.model.RT, fold.model.hpmin
    seq = decode(fold.model.al, fold.seq)
    Q = zero(T)
    p = zeros(T, n, n)
    for pt in allstruct(seq; hpmin)
        weight = exp(- energy(fold, pt) / RT)
        Q += weight
        for i = 1:n
            if isbpopening(pt, i)
                j = pt.pairs[i]
                p[i,j] += weight
                p[j,i] += weight
            end
        end
    end
    return Q, p ./ Q
end

function exhaustive_bpp_partfn(fold::Fold)
    Q, p = exhaustive_bpp_partfn(LogSR{Float64}, fold)
    logQ = Q.val
    return - fold.model.RT * logQ, map(x -> exp(x.val), p)
end


# exhaustive sequence design / inverse folding
# find sequences with highest probability of target structure
function exhaustive_design(target::Pairtable, model; nbest::Integer=20)
    if ! isvalid(target; hpmin=model.hpmin)
        throw(ArgumentError("target structure has hairpins with size < hpmin"))
    end
    best = FixedsizePQ{String,Float64}(nbest)
    for seq in allseq(target)
        ptarget = prob_of_struct(Fold(seq, model), target)
        enqueue!(best, seq, ptarget)
    end
    return collect(best)
end

exhaustive_design(dbn::AbstractString, model; nbest::Integer=20) =
    exhaustive_design(Pairtable(dbn), model; nbest)
