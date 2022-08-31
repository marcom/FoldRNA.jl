# score, energy for Fold{LoopModel}

using Printf: @printf

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

Base.@kwdef struct Extloop
    stems :: Vector{Basepair} = Basepair[]
end

Base.@kwdef struct Multiloop
    bp :: Basepair
    stems :: Vector{Basepair}
end

Base.@kwdef struct Loop
    # extloop is Basepair(0,0)
    bp :: Basepair
    # length 0: hairpin, length 1: intloop, length >= 2: multiloop
    stems :: Vector{Basepair} = Basepair[]
end

Base.:(==)(a::Loop, b::Loop) = (a.bp == b.bp && a.stems == b.stems)

function Base.show(io::IO, ::MIME"text/plain", loop::Loop)
    bp, stems = loop.bp, loop.stems
    str_stems = join(["($(bp.i), $(bp.j))" for bp in stems], ", ")
    if bp != Basepair(0,0)
        @printf(io, "    (%3d,%3d): %s\n", bp.i, bp.j, str_stems)
    else
        println(io, "      extloop: $str_stems")
    end
end

Base.@kwdef struct LoopStructure
    loops :: Vector{Loop} = [Loop(; bp=Basepair(0,0))]
end

function Base.show(io::IO, mime::MIME"text/plain", ls::LoopStructure)
    nloop = length(ls.loops)
    println(io, "LoopStructure with $nloop loop$(nloop == 1 ? "" : "s"):")
    for loop in ls.loops
        show(io, mime, loop)
    end
end

function LoopStructure(pt::Pairtable)
    function addloops_rec!(ls::LoopStructure, loop::Loop)
        for bp in loop.stems
            stems = findstems(pt, bp.i+1, bp.j-1)
            loop = Loop(bp, stems)
            push!(ls.loops, loop)
            addloops_rec!(ls, loop)
        end
    end
    n = length(pt)
    ls = LoopStructure()
    # extloop
    stems = findstems(pt, 1, n)
    append!(ls.loops[1].stems, stems)
    addloops_rec!(ls, ls.loops[1])
    return ls
end

LoopStructure(dbn::AbstractString) = LoopStructure(Pairtable(dbn))

energy(fold::Fold{M}, hairpin::Hairpin) where {M <: LoopModel} =
    score(fold, hairpin) * fold.model.unit

energy(fold::Fold{M}, intloop::Intloop) where {M <: LoopModel} =
    score(fold, intloop) * fold.model.unit

energy(fold::Fold{M}, extloop::Extloop) where {M <: LoopModel} =
    score(fold, extloop) * fold.model.unit

energy(fold::Fold{M}, multiloop::Multiloop) where {M <: LoopModel} =
    score(fold, multiloop) * fold.model.unit

function energy(fold::Fold{M}, loop::Loop; verbose::Bool=false) where {M <: LoopModel}
    bp, stems = loop.bp, loop.stems
    if bp.i == 0 && bp.j == 0
        # extloop
        en = energy(fold, Extloop(stems))
    elseif length(stems) == 0
        # hairpin
        en = energy(fold, Hairpin(bp.i, bp.j))
    elseif length(stems) == 1
        # intloop
        en = energy(fold, Intloop(bp.i, bp.j, stems[1].i, stems[1].j))
    else
        # multiloop
        en = energy(fold, Multiloop(bp, stems))
    end
    if verbose
        iobuf = IOBuffer()
        show(iobuf, MIME"text/plain"(), loop)
        str_loop = chomp(String(take!(iobuf)))
        @printf("%-30s    en = %s\n", str_loop, en)
    end
    return en
end

energy(fold::Fold{M}, ls::LoopStructure; verbose::Bool=false) where {M <: LoopModel} =
    sum(energy(fold, l; verbose) for l in ls.loops; init=zero(fold.model.unit))

energy(fold::Fold{M}, pt::Pairtable; verbose::Bool=false) where {M <: LoopModel} =
    energy(fold, LoopStructure(pt); verbose)

energy(fold::Fold{M}, dbn::String; verbose::Bool=false) where {M <: LoopModel} =
    energy(fold, Pairtable(dbn); verbose)

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
            s += m.hairpin_init[m.maxloop] + trunc(Int, m.lxc * log(len / m.maxloop))
        end
        # mismatch interior of closing basepair
        # TODO: hardcoded 3, make it part of LoopModel ?
        if len == 3
            s += score_terminal_nonGC(fold, i, j)
        elseif len > 3
            s += score_mismatch(fold, i, j, i+1, j-1, m.mismatch_hairpin)
        end
    end
    return s
end

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
            s = m.bulge_init[m.maxloop] + trunc(Int, m.lxc * log(len / m.maxloop))
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
                s = m.intloop_init[m.maxloop] + trunc(Int, m.lxc * log(len / m.maxloop))
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
                s = m.intloop_init[m.maxloop] + trunc(Int, m.lxc * log(len / m.maxloop))
            end
            s += min(m.ninio_max, abs(n1-n2) * m.ninio_m)
            s += score_mismatch_intloop(fold, i, j, k, l, b1, b2, b3, b4)
        end
    end
    return s
end

function score(fold::Fold{M}, ml::Multiloop) where {T, M <: LoopModel{T}}
    # linear multiloop energy
    m = fold.model
    nstems = length(ml.stems)
    s = zero(T)
    n = length(fold)
    # contributions for each stem
    for bp in ml.stems
        s += score_multiloop_stem(fold, bp.i, bp.j)
    end
    # closing base pair of multiloop
    i, j = ml.bp.i, ml.bp.j
    s += score_multiloop_closing_bp(fold, i, j)

    # TODO: this calculation is wrong, need to consider span of stems
    nunpaired = (j - i - 1) - 2 * nstems
    s += score_multiloop_unpaired(fold, nunpaired)

    # # Note: (nstems + 1) for the closing base pair of the multiloop
    # s += ( score_multiloop_init(fold)
    #        + nunpaired * m.multiloop_unpaired
    #        + (nstems + 1) * m.multiloop_branch )
    return s
end

# TODO: score_extloop_unpaired for unpaired bases in extloop
function score(fold::Fold{M}, extloop::Extloop) where {T, M <: LoopModel{T}}
    n = length(fold)
    s = zero(T)
    for bp in extloop.stems
        i, j = bp.i, bp.j
        s += score_extloop_stem(fold, i, j)
    end
    return s
end

# helper functions

score_terminal_nonGC(fold::Fold{M}, i::Integer, j::Integer) where {M <: LoopModel} =
    fold.model.terminal_nonGC_bp[bptype(fold, i, j)]

score_mismatch(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
               mismatch_param) where {M <: LoopModel} =
    mismatch_param[bptype(fold, i, j), fold.seq[k], fold.seq[l]]

# Note: order for l,k swapped for second base pair in some of the following functions
# TODO: is this anywhere except score_multiloop_closing_bp ? if not
#       remove above comment as there is a comment above that function

score_stack(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer) where {M <: LoopModel} =
    fold.model.stack[bptype(fold, i, j), bptype(fold, l, k)]

score_intloop11(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer) where {M <: LoopModel} =
    fold.model.intloop11[bptype(fold, i, j), bptype(fold, l, k), fold.seq[i+1], fold.seq[j-1]]

score_intloop12(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
                b1::Integer, b3::Integer, b4::Integer) where {M <: LoopModel} =
                    fold.model.intloop12[bptype(fold, i, j), bptype(fold, l, k),
                                         fold.seq[b1], fold.seq[b3], fold.seq[b4]]

score_intloop21(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
                b1::Integer, b2::Integer, b3::Integer) where {M <: LoopModel} =
                    fold.model.intloop12[bptype(fold, l, k), bptype(fold, i, j),
                                         fold.seq[b3], fold.seq[b1], fold.seq[b2]]

score_intloop22(fold::Fold{M}, i::Integer, j::Integer, k::Integer, l::Integer,
                b1::Integer, b2::Integer, b3::Integer, b4::Integer) where {M <: LoopModel} =
                    fold.model.intloop22[bptype(fold, i, j), bptype(fold, l, k),
                                         fold.seq[b1], fold.seq[b2], fold.seq[b3], fold.seq[b4]]

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

# TODO: overlap with score_multiloop_stem, simplify
function score_extloop_stem(fold::Fold{M}, i::Integer, j::Integer) where {T, M <: LoopModel{T}}
    n = length(fold)
    dangle5 = i > 1 ? i-1 : -1
    dangle3 = j < n ? j+1 : -1
    return score_stem_extloop_multiloop_dangles(fold, i, j, dangle5, dangle3,
                                                fold.model.mismatch_extloop)
end

score_extloop_unpaired(fold::Fold{M}, nunpaired::Integer) where {M <: LoopModel} =
    nunpaired * fold.model.extloop_unpaired

# TODO: remove
# score_multiloop_init(fold::Fold{M}) where {T, M <: LoopModel{T}} =
#     fold.model.multiloop_init

function score_multiloop_stem(fold::Fold{M}, i::Integer, j::Integer) where {T, M <: LoopModel{T}}
    n = length(fold)
    dangle5 = i > 1 ? i-1 : -1
    dangle3 = j < n ? j+1 : -1
    s = fold.model.multiloop_branch
    s += score_stem_extloop_multiloop_dangles(fold, i, j, dangle5, dangle3,
                                              fold.model.mismatch_multiloop)
    return s
end

function score_multiloop_closing_bp(fold::Fold{M}, i::Integer, j::Integer) where {T, M <: LoopModel{T}}
    dangle5 = i + 1
    dangle3 = j - 1
    s = fold.model.multiloop_init
    s += fold.model.multiloop_branch
    # note: reversed order of j, i, dangle3, dangle5
    s += score_stem_extloop_multiloop_dangles(fold, j, i, dangle3, dangle5,
                                              fold.model.mismatch_multiloop)
    return s
end

function score_stem_extloop_multiloop_dangles(fold::Fold{M}, i::Integer, j::Integer,
                                              dangle5::Integer, dangle3::Integer,
                                              mismatch) where {M <: LoopModel}
    # stem contributions for extloop and multiloop
    n = length(fold)
    s = score_terminal_nonGC(fold, i, j)
    if dangle5 != -1 && dangle3 != -1
        s += score_mismatch(fold, i, j, dangle5, dangle3, mismatch)
    else
        if dangle5 != -1
            s += fold.model.dangle5[bptype(fold, i, j), fold.seq[dangle5]]
        end
        if dangle3 != -1
            s += fold.model.dangle3[bptype(fold, i, j), fold.seq[dangle3]]
        end
    end
    return s
end

score_multiloop_unpaired(fold::Fold{M}, nunpaired::Integer) where {M <: LoopModel} =
    nunpaired * fold.model.multiloop_unpaired

