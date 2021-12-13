using Test
using NucleicAcidFold: Fold, Pairtable, energy, mfe, partfn, prob_of_struct,
    RNA_BPMODEL
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
    end
    @testset "LoopModel" begin
        seq = "GGGAAACCC"
        dbn = "(((...)))"
        model = LoopModel{Float64,Int,4,6,30}(alphabet=Alphabet("ACGU"))
        model.bptype .= 1
        fold = Fold(seq, model)
        @test bptype(fold, 1, 9) isa Int
        @test energy(fold, dbn) isa Quantity
    end
end
