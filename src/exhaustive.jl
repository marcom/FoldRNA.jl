function exhaustive_mfe(fold::Fold)
    n = length(fold)
    pt_min = Pairtable("."^n)
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
    # TODO: take temperature from fold.model
    RT = fold.model.RT
    Q = zero(T)
    seq = decode(fold.model.al, fold.seq)
    for pt in allstruct(seq; fold.model.hpmin)
        en = energy(fold, pt)
        Q += exp(- en / RT)
    end
    return Q
end

function exhaustive_partfn(fold::Fold)
    s = exhaustive_partfn(LogSR{Float64}, fold)
    logQ = s.val
    return - fold.model.RT * logQ
end
