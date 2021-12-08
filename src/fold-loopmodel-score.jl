# score, energy for Fold{LoopModel}

struct Hairpinloop
    i :: Int
    j :: Int
end

function isspecialhairpin(fold::Fold{M}, hairpin::Hairpinloop) where {M <: LoopModel}
    shp = fold.model.specialhairpins
    i, j = hairpin.i, hairpin.j
    return haskey(shp, @view fold.seq[i:j])
end

function score_hairpin_special(fold::Fold{M}, hairpin::Hairpinloop) where {M <: LoopModel}
    shp = fold.model.specialhairpins
    i, j = hairpin.i, hairpin.j
    return shp[@view fold.seq[i:j]]
end

function score_hairpin_init(fold::Fold{M}, len::Integer) where {M <: LoopModel}
    m = fold.model
    if len <= m.maxloop
        # small loops have precomputed values
        return m.hairpin_init[len]
    else
        # large loops get extrapolated
        return m.hairpin_init[m.maxloop] + round(Int, m.lxc * log(len / m.maxloop))
    end
end

function score_hairpin_mismatch(fold::Fold{M}, hairpin::Hairpinloop,
                                len::Integer) where {T, M <: LoopModel{T}}
    # mismatch penalty for first unpaired bases in hairpin
    # TODO
    # - hardcoded length 3 for terminal_nonGC vs mismatch penalty,
    #   take this from fold.model
    i, j = hairpin.i, hairpin.j
    if len == 3
        return score_terminal_nonGC(fold, i, j)
    elseif len > 3
        return score_mismatch(fold, i, j, i+1, j-1, fold.model.mismatch_hairpin)
    end
    return zero(T)
end

score_terminal_nonGC(fold::Fold{M}, i::Integer, j::Integer) where {M <: LoopModel} =
    fold.model.terminal_nonGC_bp[bptype(fold, i, j)]

function score_mismatch(fold::Fold{M}, i::Integer, j::Integer,
                        k::Integer, l::Integer,
                        mismatch_param) where {M <: LoopModel}
    return mismatch_param[bptype(fold, i, j), fold.seq[k], fold.seq[l]]
end

function energy(fold::Fold{M}, hairpin::Hairpinloop) where {M <: LoopModel}
    i, j = hairpin.i, hairpin.j
    len = j - i - 1
    if isspecialhairpin(fold, hairpin)
        # special hairpin loop (Triloop, Tetraloop, etc.)
        en = score_hairpin_special(fold, hairpin)
    else
        en = score_hairpin_init(fold, len) +
            score_hairpin_mismatch(fold, hairpin, len)
    end
    return en * fold.model.unit
end
