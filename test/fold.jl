using Test
using NucleicAcidFold: Fold, Pairtable, energy, mfe, partfn, prob_of_struct,
    RNA_BPMODEL, Loop, LoopStructure
using Unitful: Quantity

@testset "Fold" begin
    @testset "BpModel" begin
        seq = "GGGAAACCC"
        dbn = "(((...)))"
        pt = Pairtable(dbn)
        model = RNA_BPMODEL

        fold = Fold(seq, model)
        @test length(fold) == length(seq)
        @test energy(fold, dbn) isa Quantity
        @test energy(fold, pt) isa Quantity
        @test mfe(fold) == -9.0u"kcal/mol"
        @test partfn(fold) isa Quantity
        @test prob_of_struct(fold, dbn) isa Float64
        @test prob_of_struct(fold, pt) isa Float64
        iobuf = IOBuffer()
        show(iobuf, MIME"text/plain"(), fold)
        @test length(take!(iobuf)) > 0
    end
    @testset "LoopModel" begin
        seq = "GGGAAACCC"
        dbn = "(((...)))"
        model = LoopModel{Float64,Int,4,6,30}(alphabet=Alphabet("ACGU"))
        model.bptype .= 1
        fold = Fold(seq, model)
        @test bptype(fold, 1, 9) isa Int
        @test energy(fold, dbn) isa Quantity
        iobuf = IOBuffer()
        show(iobuf, MIME"text/plain"(), fold)
        @test length(take!(iobuf)) > 0
    end

    @testset "Loop, LoopStructure" begin
        # Loop
        loop = Loop(bp=Basepair(1,2))
        iobuf = IOBuffer()
        show(iobuf, MIME"text/plain"(), loop)
        @test length(take!(iobuf)) > 0
        @test Loop(Basepair(1,10), [Basepair(2,9)]) == Loop(Basepair(1,10), [Basepair(2,9)])
        
        # LoopStructure
        ls = LoopStructure()
        iobuf = IOBuffer()
        show(iobuf, MIME"text/plain"(), ls)
        @test length(take!(iobuf)) > 0
        pt = Pairtable("((...)(...))")
        ls = LoopStructure(pt)
        @test length(ls.loops) == 4
        @test Loop(Basepair(0, 0),  [Basepair(1, 12)]) in ls.loops
        @test Loop(Basepair(1, 12), [Basepair(2, 6), Basepair(7, 11)]) in ls.loops
        @test Loop(Basepair(2, 6),  Basepair[]) in ls.loops 
        @test Loop(Basepair(7, 11), Basepair[]) in ls.loops

        # Note: energy(::Loop) and energy(::LoopStructure) are
        # implicitly tested in test/fold-loopmodel-score.jl
    end
end
