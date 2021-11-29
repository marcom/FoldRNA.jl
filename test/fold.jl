using Test
using NucleicAcidFold: Fold, RNA_BPMODEL, Pairtable, energy, partfn, prob_of_struct
using Unitful: Quantity

@testset "Fold" begin
    seq = "GGGAAACCC"
    dbn = "(((...)))"
    pt = Pairtable(dbn)
    fold = Fold(seq, RNA_BPMODEL)
    @test energy(fold, dbn) isa Quantity
    @test energy(fold, pt) isa Quantity
    @test partfn(fold) isa Quantity
    @test prob_of_struct(fold, dbn) isa Float64
    @test prob_of_struct(fold, pt) isa Float64
end
