# score, energy for Fold{LoopModel}

struct Hairpin
    i :: Int
    j :: Int
end

score_terminal_nonGC(fold::Fold{M}, i::Integer, j::Integer) where {M <: LoopModel} =
    fold.model.terminal_nonGC_bp[bptype(fold, i, j)]

score_mismatch(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
               mismatch_param) where {M <: LoopModel} =
    mismatch_param[bptype(fold, i, j), fold.seq[k], fold.seq[l]]

function score(fold::Fold{M}, hairpin::Hairpin) where {T, M <: LoopModel{T}}
    i, j = hairpin.i, hairpin.j
    m = fold.model
    len = j - i - 1
    hpseq = @view fold.seq[i:j]
    s = zero(T)
    if haskey(m.specialhairpins, hpseq)
        # it's a special hairpin loop (Triloop, Tetraloop, etc.)
        s += m.specialhairpins[hpseq]
    else
        # hairpin init
        if len <= m.maxloop
            # small loops have precomputed values
            s += m.hairpin_init[len]
        else
            # large loops get extrapolated
            s += m.hairpin_init[m.maxloop] + round(Int, m.lxc * log(len / m.maxloop))
        end
        # mismatch interior of closing basepair
        # TODO: hardcoded 3, make it part of LoopModel
        if len == 3
            s += score_terminal_nonGC(fold, i, j)
        elseif len > 3
            s += score_mismatch(fold, i, j, i+1, j-1, m.mismatch_hairpin)
        end
    end
    return s
end

energy(fold::Fold{M}, hairpin::Hairpin) where {M <: LoopModel} =
    score(fold, hairpin) * fold.model.unit
