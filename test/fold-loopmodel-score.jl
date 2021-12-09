using Test
using NucleicAcidFold: Hairpin
using Unitful: unit

@testset "score LoopModel" begin
    @testset "hairpin" begin
        alphabet = Alphabet("ACGU")
        T = Float64
        Tseq = Int
        nb = length(alphabet)
        nbp = 6
        maxloop = 5

        # test energy, score, specialhairpins
        seq = "GGAAAACC"
        hairpin = Hairpin(2, 7)
        model = LoopModel{T,Tseq,nb,nbp,maxloop}(; alphabet)
        model.specialhairpins[encode(model, "GAAAAC")] = 42.0
        fold = Fold(seq, model)
        @test score(fold, hairpin) == 42.0
        en = energy(fold, hairpin)
        @test unit(en) == unit(model.unit)
        @test en == 42.0 * model.unit

        # test with different types and hairpin lengths
        for T in (Float64, Int)
            model = LoopModel{T,Tseq,nb,nbp,maxloop}(; alphabet)
            model.bptype = ones(Int, nb, nb)
            for hplen = 0:maxloop+10
                seq = "GG" * "A"^hplen * "CC"
                fold = Fold(seq, model)
                hairpin = Hairpin(2, 2 + hplen + 1)
                @test score(fold, hairpin) isa T
            end
        end
    end
end
