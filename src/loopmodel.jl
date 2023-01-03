using StaticArrays: MArray, @MArray
using OffsetArrays: OffsetArray

export LoopModel

"""
    LoopModel

Nearest-neighbour model (Turner model) based on scoring loops closed
by basepairs.

The implementation currently doesn't support:
  - pseudoknots
  - coaxial stacking
  - multiple strands (although there is the scalar duplex_init)
"""
Base.@kwdef mutable struct LoopModel{T,Tseq,NB,NBP,MAXLOOP}
    #          NNDB                   here
    # ---------------------   -------------------
    # WC helices              stack
    # GU pairs                stack
    # dangling ends           dangle5, dangle3
    # terminal mismatches     mismatch_*
    # hairpins                hairpin
    # bulge loops             bulge
    # internal loops          intloop*, ninio_*
    # coaxial stacking
    # multiloops              multiloop_*
    # exterior loop           extloop_unpaired
    name               :: String = ""
    alphabet           :: Alphabet
    bptype             :: MArray{Tuple{NB, NB}, Int} = @MArray zeros(Int, NB, NB)
    hpmin              :: Int = 3
    maxloop            :: Int = MAXLOOP
    RT                 :: Quantity = RT37
    RT_ustrip          :: typeof(ustrip(RT37)) = ustrip(RT37)
    # energy unit of parameters
    unit :: Quantity = 1.0u"kcal/mol"

    stack              :: MArray{Tuple{NBP, NBP}, T} = @MArray zeros(T, NBP, NBP)
    hairpin_init       :: OffsetArray{T} = OffsetArray(zeros(T, MAXLOOP+1), 0:MAXLOOP)
    bulge_init         :: OffsetArray{T} = OffsetArray(zeros(T, MAXLOOP+1), 0:MAXLOOP)
    intloop_init       :: OffsetArray{T} = OffsetArray(zeros(T, MAXLOOP+1), 0:MAXLOOP)
    intloop11          :: MArray{Tuple{NBP,NBP,NB,NB}, T} = @MArray zeros(T, NBP, NBP, NB, NB)
    intloop12          :: Array{T,5} = zeros(T, NBP, NBP, NB, NB, NB)
    intloop22          :: Array{T,6} = zeros(T, NBP, NBP, NB, NB, NB, NB)
    dangle5            :: MArray{Tuple{NBP, NB}, T} = @MArray zeros(T, NBP, NB)
    dangle3            :: MArray{Tuple{NBP, NB}, T} = @MArray zeros(T, NBP, NB)
    mismatch_hairpin   :: MArray{Tuple{NBP, NB, NB}, T} = @MArray zeros(T, NBP, NB, NB)
    mismatch_intloop   :: MArray{Tuple{NBP, NB, NB}, T} = @MArray zeros(T, NBP, NB, NB)
    mismatch_intloop1n :: MArray{Tuple{NBP, NB, NB}, T} = @MArray zeros(T, NBP, NB, NB)
    mismatch_intloop23 :: MArray{Tuple{NBP, NB, NB}, T} = @MArray zeros(T, NBP, NB, NB)
    mismatch_multiloop :: MArray{Tuple{NBP, NB, NB}, T} = @MArray zeros(T, NBP, NB, NB)
    mismatch_extloop   :: MArray{Tuple{NBP, NB, NB}, T} = @MArray zeros(T, NBP, NB, NB)
    multiloop_branch   :: T = zero(T)
    multiloop_unpaired :: T = zero(T)
    multiloop_init     :: T = zero(T)
    extloop_unpaired   :: T = zero(T)
    specialhairpins    :: Dict{Vector{Tseq}, T} = Dict{Vector{Tseq}, T}()
    ninio_m            :: T = zero(T)
    ninio_max          :: T = zero(T)
    duplex_init        :: T = zero(T)
    terminal_nonGC     :: T = zero(T)
    terminal_nonGC_bp  :: MArray{Tuple{NBP}, T} = @MArray zeros(T, NBP)
    lxc                :: Float64 = zero(Float64)
end

encode(m::LoopModel, iter) = encode(m.alphabet, iter)
decode(m::LoopModel, iter) = decode(m.alphabet, iter)
bptype(m::LoopModel, ci::Integer, cj::Integer) = m.bptype[ci, cj]
canbp(m::LoopModel, ci::Integer, cj::Integer) = bptype(m, ci, cj) != 0

"""
    filter_wildcard_chars(model::LoopModel)

Remove all parameters from tables that refer to wildcard nucleotides,
e.g. 'N'.
"""
function filter_wildcard_chars(oldm::LoopModel{T,Tseq,NB_old,NBP_old,MAXLOOP}) where {T,Tseq,NB_old,NBP_old,MAXLOOP}
    # TODO: use nospecialize on LoopModel arg?
    # copy that also works for immutables
    # see: https://github.com/JuliaLang/julia/issues/45143
    function mycopy(x::T) where T
        if ismutable(x) && ! (x isa String)
            return copy(x)
        else
            return x
        end
    end

    # TODO: assert length(old_chars) == NB_old
    # assert maximum(oldm.bptype) == NBP
    old_chars = oldm.alphabet.chars
    old_wc = oldm.alphabet.wildcard_chars
    idxs_wc = [i for i = 1:NB_old if old_chars[i] in old_wc]
    new_chars = [c for c in old_chars if !(c in old_wc)]
    new_wc = Char[]

    # map old indices to new indices (in alphabet)
    map_idxs = zeros(Int, length(old_chars))
    for (i, c) in enumerate(old_chars)
        k = findfirst(x -> x == c, new_chars)
        map_idxs[i] = isnothing(k) ? 0 : k
    end
    Ib = [i for i in 1:NB_old if !(old_chars[i] in old_wc)]

    old_basepairs = Tuple{Int,Int}[]
    for i = 1:NBP_old
        idxs = findfirst(x -> x == i, oldm.bptype)
        !isnothing(idxs) && push!(old_basepairs, idxs)
    end
    new_basepairs = [(map_idxs[i],map_idxs[j]) for (i,j) in old_basepairs if !(i in idxs_wc || j in idxs_wc)]
    Ibp = Int[i for i in 1:min(NBP_old,length(old_basepairs)) if old_basepairs[i] in new_basepairs]

    NB = length(new_chars)
    NBP = length(new_basepairs)
    alphabet = Alphabet(oldm.alphabet.name, new_chars; wildcard_chars=new_wc)
    m = LoopModel{T,Tseq,NB,NBP,MAXLOOP}(; alphabet)

    # copy over everything that doesn't have to be filtered
    fields_to_be_filtered = [
        :alphabet, :bptype,
        :stack, :intloop11, :intloop12, :intloop22, :dangle5, :dangle3,
        :mismatch_hairpin, :mismatch_intloop, :mismatch_intloop1n,
        :mismatch_intloop23, :mismatch_multiloop, :mismatch_extloop,
        :terminal_nonGC_bp
    ]
    for field in setdiff(fieldnames(typeof(oldm)), fields_to_be_filtered)
        setfield!(m, field, mycopy(getfield(oldm, field)))
    end

    # filter properties with potential wildcard chars
    for (k, (i,j)) in enumerate(new_basepairs)
        m.bptype[i,j] = k
    end
    m.stack .= @view oldm.stack[Ibp, Ibp]
    m.intloop11 .= @view oldm.intloop11[Ibp, Ibp, Ib, Ib]
    m.intloop12 .= @view oldm.intloop12[Ibp, Ibp, Ib, Ib, Ib]
    m.intloop22 .= @view oldm.intloop22[Ibp, Ibp, Ib, Ib, Ib, Ib]
    m.dangle5 .= @view oldm.dangle5[Ibp, Ib]
    m.dangle3 .= @view oldm.dangle3[Ibp, Ib]
    m.mismatch_hairpin .= @view oldm.mismatch_hairpin[Ibp, Ib, Ib]
    m.mismatch_intloop .= @view oldm.mismatch_intloop[Ibp, Ib, Ib]
    m.mismatch_intloop1n .= @view oldm.mismatch_intloop1n[Ibp, Ib, Ib]
    m.mismatch_intloop23 .= @view oldm.mismatch_intloop23[Ibp, Ib, Ib]
    m.mismatch_multiloop .= @view oldm.mismatch_multiloop[Ibp, Ib, Ib]
    m.mismatch_extloop .= @view oldm.mismatch_extloop[Ibp, Ib, Ib]
    m.terminal_nonGC_bp .= @view oldm.terminal_nonGC_bp[Ibp]

    return m
end


# loopmodel CKY algorithm
# T must be a semiring (TODO: list exact properties used)
function loopmodel(T::Type, to_SR::Function, fold; hpmin::Integer=3, maxintloop::Integer=-1)
    h = hpmin
    n = length(fold)
    # TODO: use upper-diagonal matrices to save space
    A   = zeros(T, n, n)  # property A on subsequence [i,j] (extloop)
    Ab  = zeros(T, n, n)  # like A, but (i,j) form a base pair
    Am  = zeros(T, n, n)  # like A, but part of a multiloop with at least one stem
    Am1 = zeros(T, n, n)  # like A, but part of a multiloop with exactly one stem from basepair (i,k), i<k<=j

    # base case init for A (extloop)
    for d = 0:h
        len = d + 1
        for i = 1:n-d
            j = i + d
            A[i,j] = to_SR(score_extloop_unpaired(fold, len))

            # TODO: remove
            #A[i,j] = one(T)
        end
    end

    for d = h+1:n-1
        for i = 1:n-d
            j = i + d
            # Ab: (i,j) form a base pair
            Ab[i,j] = loopmodel_Ab(to_SR, fold, i, j, hpmin, maxintloop, Ab, Am, Am1)
            # Am: part of a multiloop with at least one stem
            Am[i,j] = loopmodel_Am(to_SR, fold, i, j, hpmin, maxintloop, Am, Ab)
            # Am1: part of a multiloop with exactly one stem from basepair (i,k), i<k<=j
            Am1[i,j] = loopmodel_Am1(to_SR, fold, i, j, hpmin, maxintloop, Am1, Ab)
            # A: extloop
            A[i,j] = loopmodel_A(to_SR, fold, i, j, hpmin, maxintloop, A, Ab)
        end
    end
    return A, Ab, Am, Am1 # A[1,n]
end

function loopmodel_Ab(to_SR, fold, i, j, hpmin, maxintloop, Ab, Am, Am1)
    canbp(fold, i, j) || return zero(eltype(Ab))
    h = hpmin
    # Ab: (i,j) form a base pair
    # case: (i,j) closes hairpin
    ab_ij = to_SR(score(fold, Hairpin(i, j)))
    # case: (i,j) closes internal loop started from (k,l)
    if maxintloop >= 0
        for k = i+1:min(j-h-2, i + maxintloop + 1)
            # k-i-1: unpaired bases before k
            # j-i-1 - (l-k+1): unpaired bases in loop
            #    = j - (l-i) + k - 2
            #for l = k+h+1:j-1
            for l = max(k+h+1, k + j - i - 2 - maxintloop):j-1
                loopsize = j-i-1 - (l-k+1)
                # maxintloop = j-i-1 - (l-k+1)
                # <=> maxintloop - j + i + 1 = -l + k - 1
                # <=> l = k + j - i - 2 - maxintloop
                if loopsize <= maxintloop
                    if canbp(fold, k, l)
                        ab_ij += to_SR(score(fold, Intloop(i, j, k, l))) * Ab[k,l]
                    end
                else
                    println("in loopmodel, intloop: unnecessary iteration: i = $i, j= $j, k = $k, l = $l")
                end
            end
        end
    else
        for k = i+1:j-h-2
            for l = k+h+1:j-1
                if canbp(fold, k, l)
                    ab_ij += to_SR(score(fold, Intloop(i, j, k, l))) * Ab[k,l]
                end
            end
        end
    end
    # case: (i,j) closes multiloop
    for k = i+h+2:j-h-3
        ab_ij += Am[i+1,k] * Am1[k+1,j-1] * to_SR(score_multiloop_closing_bp(fold, i, j))

        # TODO: remove
        # ab_ij += Am[i+1,k] * Am1[k+1,j-1] *
        #     to_SR(score_multiloop_closing_bp(fold, i, j)) *
        #     to_SR(score_multiloop_init(fold))
    end
    return ab_ij
end

function loopmodel_Am(to_SR, fold, i, j, hpmin, maxintloop, Am, Ab)
    h = hpmin
    # Am: part of a multiloop with at least one stem
    # case: j unpaired
    am_ij = Am[i,j-1] * to_SR(score_multiloop_unpaired(fold, 1))
    # case: j paired to k (i <= k < j-h), exactly one stem in multiloop
    for k = i:j-h-1
        # TODO: unnecessary check for canbp if Ab[k,j] == zero(eltype(Ab)) when (k,j) can't basepair
        if canbp(fold, k, j)
            am_ij += to_SR(score_multiloop_unpaired(fold, k - i)) * Ab[k,j] * to_SR(score_multiloop_stem(fold, k, j))
        end
    end
    # case: j paired to k (i+h+2 <= k < j-h), more than one stem in multiloop
    for k = i+h+2:j-h-1
        # TODO: unnecessary check for canbp if Ab[k,j] == zero(eltype(Ab)) when (k,j) can't basepair
        if canbp(fold, k, j)
            am_ij += Am[i,k-1] * Ab[k,j] * to_SR(score_multiloop_stem(fold, k, j))
        end
    end
    return am_ij
end

function loopmodel_Am1(to_SR, fold, i, j, hpmin, maxintloop, Am1, Ab)
    h = hpmin
    # Am1: part of a multiloop with exactly one stem from basepair (i,k), i<k<=j
    # case: j unpaired
    am1_ij = Am1[i,j-1] * to_SR(score_multiloop_unpaired(fold, 1))
    # case: j paired to i (must be to i due to definition of Am1)
    # TODO: unnecessary check for canbp if Ab[i,j] == zero(eltype(Ab)) when (i,j) can't basepair
    if canbp(fold, i, j)
        am1_ij += Ab[i,j] * to_SR(score_multiloop_stem(fold, i, j))
    end
    return am1_ij
end

function loopmodel_A(to_SR, fold, i, j, hpmin, maxintloop, A, Ab)
    h = hpmin
    # A: extloop

    # TODO: case distinction done on i here, but j elsewhere
    #       was there a good reason for that? or just arbitrary?
    # case: i unpaired
    a_ij = to_SR(score_extloop_unpaired(fold, 1)) * A[i+1,j]
    # case: (i,k) paired, i+h+1 <= k < j
    for k = i+h+1:j-1
        # TODO: unnecessary check for canbp if Ab[i,k] == zero(eltype(Ab)) when (i,k) can't basepair
        if canbp(fold, i, k)
            a_ij += Ab[i,k] * to_SR(score_extloop_stem(fold, i, k)) * A[k+1,j]
        end
    end
    # case: (i,k) paired, k == j
    # TODO: unnecessary check for canbp if Ab[i,j] == zero(eltype(Ab)) when (i,j) can't basepair
    if canbp(fold, i, j)
        a_ij += Ab[i,j] * to_SR(score_extloop_stem(fold, i, j))
    end
    return a_ij
end
