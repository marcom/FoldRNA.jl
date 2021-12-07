using Test
using NucleicAcidFold: Fold, RNA_BPMODEL, Pairtable, energy, mfe, partfn,
    prob_of_struct
using Unitful: Quantity

@testset "Fold" begin
    @testset "for BpModel" begin
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
    @testset "for LoopModel" begin
        seq = "GGGAAACCC"
        model = LoopModel{Float64,Int,4,6,30}(alphabet=Alphabet("ACGU"))
        fold = Fold(seq, model)
        @test bptype(fold, 1, 9) isa Int
    end
end
