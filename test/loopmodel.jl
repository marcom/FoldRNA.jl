using Test
using NucleicAcidFold

@testset "loopmodel" begin
    @testset "LoopModel" begin
        @test LoopModel{Float64,Int,4,6,30}(alphabet=Alphabet("ACGU")) isa LoopModel
    end
end
