using Test
using FoldRNA: bptype, canbp
using Unitful: Quantity

@testset "loopmodel" begin
    @testset "LoopModel" begin
        model = LoopModel{Float64,Int,4,6,30}(alphabet=Alphabet("ACGU"))
        @test model isa LoopModel
        @test bptype(model, 1, 1) isa Int
        seq = "GGGAAACCC"
        fc = Fold(seq, model)
        @test mfe(fc) isa Quantity
        @test partfn(fc) isa Quantity

        model = RNA_TURNER2004
        @test model isa LoopModel
        @test bptype(model, 1, 1) isa Int
        @test canbp(model, 1, 4)
        @test ! canbp(model, 1, 2)
        seq = "GGGAAACCC"
        fc = Fold(seq, model)
        @test mfe(fc) isa Quantity
        @test partfn(fc) isa Quantity
    end
end
