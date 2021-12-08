using Test
using NucleicAcidFold: Hairpinloop, isspecialhairpin, score_hairpin_special,
    score_hairpin_init, score_hairpin_mismatch
using Unitful: unit

@testset "score LoopModel" begin
    @testset "hairpin" begin
        alphabet = Alphabet("ACGU")
        T = Float64
        Tseq = Int
        nb = length(alphabet)
        nbp = 6
        maxloop = 5

        seq = "GGAAAACC"
        hairpin = Hairpinloop(2, 7)
        model = LoopModel{T,Tseq,nb,nbp,maxloop}(; alphabet)
        model.specialhairpins[encode(model, "GAAAAC")] = 42.0
        fold = Fold(seq, model)
        @test isspecialhairpin(fold, hairpin)
        @test score_hairpin_special(fold, hairpin) == 42.0
        en = energy(fold, hairpin)
        @test unit(en) == unit(model.unit)
        @test en == 42.0 * model.unit

        for T in (Float64, Int)
            model = LoopModel{T,Tseq,nb,nbp,maxloop}(; alphabet)
            model.bptype = ones(Int, nb, nb)
            for hplen = maxloop+10:maxloop+10
                seq = "GG" * "A"^hplen * "CC"
                fold = Fold(seq, model)
                hairpin = Hairpinloop(2, 2 + hplen + 1)
                @test score_hairpin_init(fold, hplen) isa T
                @test score_hairpin_mismatch(fold, hairpin, hplen) isa T
            end
        end
    end
end
