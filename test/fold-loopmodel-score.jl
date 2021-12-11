using Test
using NucleicAcidFold: Hairpin, Intloop, Multiloop, Extloop
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
        for T in (Int, Float64)
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

    @testset "intloop" begin
        for T in (Int, Float64)
            #      1234567890123456789
            seq = "GGGGGGGGAAACCCCCCCC"
            model = LoopModel{T,Int,4,6,30}(alphabet=Alphabet("ACGU"))
            model.bptype = ones(Int, 4, 4)
            fold = Fold(seq, model)
            for i = 1:1, j = 19:19, k = i+1:8, l = 12:j-1
                intloop = Intloop(i, j, k, l)
                @test score(fold, intloop) isa T
                @test unit(energy(fold, intloop)) == unit(model.unit)
            end
        end
    end

    @testset "multiloop" begin
        TESTS_MULTILOOP = [
            # 12345678901234567
            ("AGAGAAACAGAAACACA", Multiloop(bp=Basepair(2,16), stems=[Basepair(4,8), Basepair(10,14)])),
            # .(.(...).(...).).
            ("GAAGAAACAGAAACAAC", Multiloop(bp=Basepair(1,17), stems=[Basepair(4,8), Basepair(10,14)])),
            # (..(...).(...)..)
            ("AGAGAAACAGAAACAAC", Multiloop(bp=Basepair(2,17), stems=[Basepair(4,8), Basepair(10,14)])),
            # .(.(...).(...)..)
            ("GAAGAAACAGAAACACA", Multiloop(bp=Basepair(1,16), stems=[Basepair(4,8), Basepair(10,14)])),
            # (..(...).(...).).
        ]
        for (seq, multiloop) in TESTS_MULTILOOP
            for T in (Int, Float64)
                model = LoopModel{T,Int,4,6,30}(alphabet=Alphabet("ACGU"))
                model.bptype = ones(Int, 4, 4)
                fold = Fold(seq, model)
                @test score(fold, multiloop) isa T
                @test unit(energy(fold, multiloop)) == unit(model.unit)
            end
        end
    end

    @testset "extloop" begin
        for T in (Int, Float64)
            seq = "AGAAACAGAAACG"
            # ml   .(...).(...).
            #      1234567890123
            model = LoopModel{T,Int,4,6,30}(alphabet=Alphabet("ACGU"))
            model.bptype = ones(Int, 4, 4)
            fold = Fold(seq, model)
            extloop = Extloop(stems=[Basepair(2,6), Basepair(8,12)])
            @test score(fold, extloop) isa T
            @test unit(energy(fold, extloop)) == unit(model.unit)
        end
    end
end
