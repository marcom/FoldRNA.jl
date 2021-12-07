using Test
using NucleicAcidFold: bptype

@testset "loopmodel" begin
    @testset "LoopModel" begin
        model = LoopModel{Float64,Int,4,6,30}(alphabet=Alphabet("ACGU"))
        @test model isa LoopModel
        @test bptype(model, 1, 1) isa Int
    end
end
