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
    maxloop = MAXLOOP
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
    # energy unit of parameters
    unit :: Quantity = 1.0u"kcal/mol"
end

encode(m::LoopModel, iter) = encode(m.alphabet, iter)
decode(m::LoopModel, iter) = decode(m.alphabet, iter)
bptype(m::LoopModel, si::Integer, sj::Integer) = m.bptype[si, sj]
