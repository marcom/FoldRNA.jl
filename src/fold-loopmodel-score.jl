# score, energy for Fold{LoopModel}

struct Hairpin
    i :: Int
    j :: Int
end

struct Intloop
    i :: Int
    j :: Int
    k :: Int
    l :: Int
end

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

function score(fold::Fold{M}, intloop::Intloop) where {T, M <: LoopModel{T}}
    i, j, k, l = intloop.i, intloop.j, intloop.k, intloop.l
    m = fold.model
    s = zero(T)
    n1 = k - i - 1
    n2 = j - l - 1
    b1 = i + 1
    b2 = k - 1
    b3 = l + 1
    b4 = j - 1
    if (n1, n2) == (0,0)
        # stacked base pairs
        s = score_stack(fold, i, j, k, l)
    elseif n1 == 0 || n2 == 0
        # bulge loop
        len = max(n1, n2)
        # initiation penalty
        if len <= m.maxloop
            s = m.bulge_init[len]
        else
            # large loops get extrapolated
            s = m.bulge_init[m.maxloop] + round(Int, m.lxc * log(len / m.maxloop))
        end
        if n1 == 1 || n2 == 1
            # short 0×1 bulge loop, so it gets stacking
            s += score_stack(fold, i, j, k, l)
        else
            s += score_terminal_nonGC(fold, i, j)
            s += score_terminal_nonGC(fold, k, l)
        end
    else
        # interior loop
        if (n1,n2) == (1,1)
            # 1×1 interior loop
            s = score_intloop11(fold, i, j, k, l)
        elseif (n1,n2) == (1,2)
            # 1×2 interior loop
            s = score_intloop12(fold, i, j, k, l, b1, b3, b4)
        elseif (n1,n2) == (2,1)
            # 2×1 interior loop
            s = score_intloop21(fold, i, j, k, l, b1, b2, b3)
        elseif n1 == 1 || n2 == 1
            # 1×n interior loop
            # TODO: large overlap with generic loop calc
            len = n1 + n2
            # initiation penalty
            if len <= m.maxloop
                s = m.intloop_init[len]
            else
                # large loops get extrapolated
                s = m.intloop_init[m.maxloop] + round(Int, m.lxc * log(len / m.maxloop))
            end
            s += min(m.ninio_max, abs(n1-n2) * m.ninio_m)
            s += score_mismatch_intloop1n(fold, i, j, k, l, b1, b2, b3, b4)
        elseif (n1,n2) == (2,2)
            # 2×2 interior loop
            s = score_intloop22(fold, i, j, k, l, b1, b2, b3, b4)
        elseif (n1,n2) == (2,3) || (n1,n2) == (3,2)
            # 2×3 interior loop
            # TODO: factor out ninio calc
            s = m.intloop_init[5] + min(m.ninio_max, abs(n1 - n2) * m.ninio_m)
            s += score_mismatch_intloop23(fold, i, j, k, l, b1, b2, b3, b4)
        else
            # generic interior loop
            # TODO: large overlap with 1×n loop calc
            len = n1 + n2
            # initiation penalty
            if len <= m.maxloop
                s = m.intloop_init[len]
            else
                # large loops get extrapolated
                s = m.intloop_init[m.maxloop] + round(Int, m.lxc * log(len / m.maxloop))
            end
            s += min(m.ninio_max, abs(n1-n2) * m.ninio_m)
            s += score_mismatch_intloop(fold, i, j, k, l, b1, b2, b3, b4)
        end
    end
    return s
end

energy(fold::Fold{M}, intloop::Intloop) where {T, M <: LoopModel{T}} =
    score(fold, intloop) * fold.model.unit


# helper functions

score_terminal_nonGC(fold::Fold{M}, i::Integer, j::Integer) where {M <: LoopModel} =
    fold.model.terminal_nonGC_bp[bptype(fold, i, j)]

score_mismatch(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
               mismatch_param) where {M <: LoopModel} =
    mismatch_param[bptype(fold, i, j), fold.seq[k], fold.seq[l]]

# Note: order for l,k swapped for second base pair in some of the following functions

score_stack(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer) where {M <: LoopModel} =
    fold.model.stack[bptype(fold, i, j), bptype(fold, l, k)]

score_intloop11(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer) where {M <: LoopModel} =
    fold.model.intloop11[bptype(fold, i, j), bptype(fold, l, k), fold.seq[i+1], fold.seq[j-1]]

score_intloop12(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
                b1::Integer, b3::Integer, b4::Integer) where {M <: LoopModel} =
    fold.model.intloop12[bptype(fold, i, j), bptype(fold, l, k), fold.seq[b1], fold.seq[b3], fold.seq[b4]]

score_intloop21(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
                b1::Integer, b2::Integer, b3::Integer) where {M <: LoopModel} =
    fold.model.intloop12[bptype(fold, l, k), bptype(fold, i, j), fold.seq[b3], fold.seq[b1], fold.seq[b2]]

score_intloop22(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
                b1::Integer, b2::Integer, b3::Integer, b4::Integer) where {M <: LoopModel} =
    fold.model.intloop22[bptype(fold, i, j), bptype(fold, l, k), fold.seq[b1], fold.seq[b2], fold.seq[b3], fold.seq[b4]]

score_mismatch_intloop1n(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
                         b1::Integer, b2::Integer, b3::Integer, b4::Integer) where {M <: LoopModel} =
    fold.model.mismatch_intloop1n[bptype(fold, i, j), fold.seq[b1], fold.seq[b4]] +
        fold.model.mismatch_intloop1n[bptype(fold, l, k), fold.seq[b3], fold.seq[b2]]

score_mismatch_intloop(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
                       b1::Integer, b2::Integer, b3::Integer, b4::Integer) where {M <: LoopModel} =
    fold.model.mismatch_intloop[bptype(fold, i, j), fold.seq[b1], fold.seq[b4]] +
        fold.model.mismatch_intloop[bptype(fold, l, k), fold.seq[b3], fold.seq[b2]]

score_mismatch_intloop23(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
                         b1::Integer, b2::Integer, b3::Integer, b4::Integer) where {M <: LoopModel} =
    fold.model.mismatch_intloop23[bptype(fold, i, j), fold.seq[b1], fold.seq[b4]] +
        fold.model.mismatch_intloop23[bptype(fold, l, k), fold.seq[b3], fold.seq[b2]]
