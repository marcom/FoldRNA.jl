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
bptype(m::LoopModel, si::Integer, sj::Integer) = m.bptype[si, sj]


# loopmodel CKY algorithm
function loopmodel(T::Type, fold; hpmin::Integer=3, maxintloop::Integer=-1)
    h = hpmin
    n = length(fold)
    # TODO: use upper-diagonal matrices to save space
    A   = zeros(T, n, n)  # property A on subsequence [i,j] (extloop)
    Ab  = zeros(T, n, n)  # like A, but (i,j) form a base pair
    Am  = zeros(T, n, n)  # like A, but part of a multiloop with at least one stem
    Am1 = zeros(T, n, n)  # like A, but part of a multiloop with exactly one stem with (i,h) the closing base pair (i<h<=j)

    # base case initialisation
    for d = 0:h
        for i = 1:n-d
            j = i + d
            A[i,j] = one(T)
        end
    end

    for d = h+1:n-1
        for i = 1:n-d
            j = i + d
            # Ab: (i,j) form a base pair
            if canbp(fold, i, j)
                # case: (i,j) closes hairpin
                Ab[i,j] = score(fold, Hairpin(i, j))
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
                                    Ab[i,j] += score_intloop(fold, i, j, k, l) * Ab[k,l]
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
                                Ab[i,j] += score_intloop(fold, i, j, k, l) * Ab[k,l]
                            end
                        end
                    end
                end
                # case: (i,j) closes multiloop
                for k = i+h+2:j-h-3
                    Ab[i,j] += Am[i+1,k] * Am1[k+1,j-1] * score_multiloop_closing_stem(fold, i, j)
                end
            end

            # Am: part of a multiloop with at least one stem
            # case: j unpaired
            Am[i,j] = Am[i,j-1] * score_multiloop_unpaired(fold, 1)
            # case: j paired to k (i <= k < j-h), exactly one stem in multiloop
            for k = i:j-h-1
                # TODO: unnecessary check for canbp if Ab[k,j] == 0 if (k,j) can't basepair
                if canbp(fold, k, j)
                    Am[i,j] += score_multiloop_unpaired(fold, k - i) * Ab[k,j]
                end
            end
            # case: j paired to k (i+h+2 <= k < j-h), more than one stem in multiloop
            for k = i+h+2:j-h-1
                # TODO: unnecessary check for canbp if Ab[k,j] == 0 if (k,j) can't basepair
                if canbp(fold, k, j)
                    Am[i,j] += Am[i,k-1] * Ab[k,j] * score_multiloop_stem(fold, k, j)
                end
            end

            # Am1: part of a multiloop with exactly one stem with (i,h) the closing base pair (i<h<=j)
            # case: j unpaired
            Am1[i,j] =  Am1[i,j-1] * score_multiloop_unpaired(fold, 1)
            # case: j paired to i (must be to i due to definition of Am1)
            # TODO: unnecesary check for canbp if Ab[i,j] == 0 if (i,j) can't basepair
            if canbp(fold, i, j)
                Am1[i,j] += Ab[i,j] * score_multiloop_stem(fold, i, j)
            end

            # A: extloop
            # case: i unpaired
            A[i,j] = A[i+1,j] * score_extloop_unpaired(fold, 1)
            # case: (i,k) paired, i+h+1 <= k < j
            for k = i+h+1:j-1
                # TODO: unnecesary check for canbp if Ab[i,k] == 0 if (i,k) can't basepair
                if canbp(fold, i, k)
                    A[i,j] += Ab[i,k] * A[k+1,j]
                end
            end
            # case: (i,k) paired, k == j
            # TODO: unnecesary check for canbp if Ab[i,k] == 0 if (i,k) can't basepair
            if canbp(fold, i, j)
                A[i,j] += Ab[i,j]
            end
        end
    end
    return A, Ab, Am, Am1 # A[1,n]
end

