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
