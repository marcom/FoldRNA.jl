using Test
using NucleicAcidFold: FoldPseq, RNA_BPMODEL, Pairtable, energy, partfn,
    prob_of_struct, onehot, isohot, hotidx
using Unitful: Quantity

@testset "FoldPseq" begin
    seq = "GGGAAACCC"
    dbn = "(((...)))"
    pt = Pairtable(dbn)
    fold = FoldPseq(seq, RNA_BPMODEL)
    @test length(fold) == length(seq)
    @test String(fold) == seq
    @test energy(fold, dbn) isa Quantity
    @test energy(fold, pt) isa Quantity
    @test partfn(fold) isa Quantity
    @test prob_of_struct(fold, dbn) isa Float64
    @test prob_of_struct(fold, pt) isa Float64

    @testset "onehot" begin
        T = Float16
        chars = "ABC"
        seq = "CBACCA"
        al = Alphabet(chars)

        x = onehot(T, chars, seq)
        @test x isa Matrix{T}
        @test x == [ 0. 0. 1. 0. 0. 1.
                     0. 1. 0. 0. 0. 0.
                     1. 0. 0. 1. 1. 0. ]
        @test onehot(chars, seq) isa Matrix{Float64}
        @test onehot(chars, seq) == onehot(Float64, chars, seq)

        @test onehot(al, seq) isa Matrix{Float64}
        @test onehot(al, seq) == x
        @test onehot(T, al, seq) isa Matrix{T}
        @test onehot(T, al, seq) == x
    end
    @testset "isohot" begin
        for T in (Float16, Float64, BigFloat)
            for (nd, n) in ((1,1), (2,3), (3,5))
                @test isohot(T, nd, n) isa Matrix{T}
                x = isohot(nd, n)
                @test x isa Matrix{Float64}
                @test all(i -> i == 1/nd, x)

                al = Alphabet(join('A' + i for i = 0:nd-1))
                @test isohot(T, al, n) isa Matrix{T}
                x = isohot(al, n)
                @test x isa Matrix{Float64}
                @test all(i -> i == 1/nd, x)
            end
        end
    end
    @testset "hotidx" begin
        al = Alphabet("ABC")
        seq = ["B", "CAAB", "BACBACC"]
        for s in seq
            o = onehot(al, s)
            h = hotidx(o)
            @test h isa Vector{Int}
            @test h == hotidx(o; dims=1) == encode(al, s)
        end
    end
end
