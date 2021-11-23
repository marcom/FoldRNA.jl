using Test
using Unitful: Quantity
using NucleicAcidFold: exhaustive_partfn

@testset "exhaustive" begin
    @testset "bpmodel partfn" begin
        seq = "GGGAAACCC"
        param = NucleicAcidFold.DEFAULT_BPMODEL_PARAM
        for T in (Float64, BigFloat, LogSR{Float64})
            @test exhaustive_partfn(T, seq, param) isa T
            @test exhaustive_partfn(seq, NucleicAcidFold.DEFAULT_BPMODEL_PARAM) isa Quantity
        end
    end
end
