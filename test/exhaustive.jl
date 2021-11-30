using Test
using Unitful: Quantity
using NucleicAcidFold: Fold, exhaustive_mfe, exhaustive_partfn, RNA_BPMODEL

@testset "exhaustive" begin
    @testset "mfe for bpmodel" begin
        seq = "GGGAAACCC"
        model = RNA_BPMODEL
        @test exhaustive_mfe(Fold(seq, model)) == (-9.0u"kcal/mol", Pairtable("(((...)))"))
    end
    @testset "partfn for bpmodel" begin
        model = RNA_BPMODEL
        seq = "GGGAAACCC"
        for T in (Float64, BigFloat, LogSR{Float64})
            fold = Fold(seq, model)
            @test exhaustive_partfn(T, fold) isa T
            @test exhaustive_partfn(fold) isa Quantity
        end

        for seq in ["UUUCGAAGUUAGUCA", "CGUGGUCCUCUCCGU"]
            fold = Fold(seq, model)
            RTlogQ = model.RT * log(exhaustive_partfn(Float64, fold))
            @test -RTlogQ ≈ partfn(fold)
            @test exhaustive_partfn(fold) ≈ partfn(fold)
        end
    end
end
