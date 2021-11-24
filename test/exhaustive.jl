using Test
using Unitful: Quantity
using NucleicAcidFold: exhaustive_partfn

@testset "exhaustive" begin
    @testset "bpmodel partfn" begin
        param = NucleicAcidFold.DEFAULT_BPMODEL_PARAM
        seq = "GGGAAACCC"
        for T in (Float64, BigFloat, LogSR{Float64})
            @test exhaustive_partfn(T, seq, param) isa T
            @test exhaustive_partfn(seq, param) isa Quantity
        end

        param = NucleicAcidFold.DEFAULT_BPMODEL_PARAM
        for seq in ["UUUCGAAGUUAGUCA", "CGUGGUCCUCUCCGU"]
            @test exhaustive_partfn(seq, param) ≈ partfn(seq, param)
        end
    end
end
