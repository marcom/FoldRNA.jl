using Test
using FoldRNA: Basepair, score, score_exp

@testset "score BpModel" begin
    showtestset()
    seq = "GGGAAACCC"
    model = RNA_BPMODEL
    fold = Fold(seq, model)
    @test score(fold, Basepair(2, 8)) isa Number
    @test score_exp(fold, Basepair(2, 8)) isa Number
end

