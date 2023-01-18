using Test
using FoldRNA: bptype, canbp
using Unitful: Quantity

@testset "loopmodel" begin
    showtestset()
    @testset "LoopModel" begin
        model = LoopModel{Float64,Int,4,6,30}(alphabet=Alphabet("ACGU"))
        @test model isa LoopModel
        @test bptype(model, 1, 1) isa Int
        seq = "GGGAAACCC"
        fc = Fold(seq, model)
        @test mfe(fc) isa Quantity
        @test partfn(fc) isa Quantity
        m2 = FoldRNA.filter_wildcard_chars(model)
        @test m2 isa LoopModel

        model = RNA_TURNER2004
        @test model isa LoopModel
        @test bptype(model, 1, 1) isa Int
        @test canbp(model, 1, 4)
        @test ! canbp(model, 1, 2)
        seq = "GGGAAACCC"
        fc = Fold(seq, model)
        @test mfe(fc) isa Quantity
        @test partfn(fc) isa Quantity
        m2 = FoldRNA.filter_wildcard_chars(model)
        @test m2 isa LoopModel
        @test Set(m2.alphabet.chars) == Set(['A', 'C', 'G', 'U'])
        # TODO: this depends on ('N','N') being the last basepair and
        # the only with wildcard chars
        @test all(m2.stack .== @view model.stack[1:end-1, 1:end-1])
        # TODO: this depends on 'N' being the last base and ('N','N')
        # being the last basepair
        @test all(m2.intloop11 .== @view model.intloop11[1:end-1, 1:end-1, 1:end-1, 1:end-1])

        # Base.print
        iobuf = IOBuffer()
        print(iobuf, RNA_TURNER2004)
        @test String(take!(iobuf)) == "LoopModel: rna_turner2004.par"
    end
end
