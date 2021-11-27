using Test
using NucleicAcidFold: allseq, numseq, DEFAULT_BASES, DEFAULT_BASEPAIRS

@testset "allseq" begin
    for n = 1:5
        @test sum(1 for _ in allseq(n)) == numseq("."^n)
    end
    for (bases, basepairs) in [
        DEFAULT_BASES => DEFAULT_BASEPAIRS,
        "A"   => ["AA"],
        "AB"  => ["AA", "BB"],
        "ABC" => ["AB", "BA"],
        ]
        nb = length(bases)
        nbp = length(basepairs)
        for n = 1:6
            @test sum(1 for _ in allseq(n; bases)) == numseq("."^n; nbases=nb)
        end
        for dbn in ["...", "(.)", ".((..))"]
            nseq = numseq(dbn; nbases=nb, nbasepairs=nbp)
            @test sum(1 for _ in allseq(dbn; bases, basepairs)) == nseq
            @test length(unique(collect(allseq(dbn; bases, basepairs)))) == nseq
        end
    end
end
