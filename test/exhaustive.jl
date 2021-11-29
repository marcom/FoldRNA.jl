using Test
using Unitful: Quantity
using NucleicAcidFold: Fold, exhaustive_partfn, RNA_BPMODEL

@testset "exhaustive" begin
    @testset "bpmodel partfn" begin
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
