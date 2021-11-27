using Test
using NucleicAcidFold: allseq, numseq, DEFAULT_BASES, DEFAULT_BASEPAIRS

const TEST_BASE_BASEPAIRS = [
    DEFAULT_BASES => DEFAULT_BASEPAIRS,
    "A"   => ["AA"],
    "AB"  => ["AA", "BB"],
    "ABC" => ["AB", "BA"],
]

const TEST_STRUCTURES_DBN = [
    "...",
    "(.)",
    ".((..))",
]


@testset "allseq" begin
    @testset "allseq(n)" begin
        for n = 1:5
            @test sum(1 for _ in allseq(n)) == numseq("."^n)
        end
        for (bases, basepairs) in TEST_BASE_BASEPAIRS
            nb = length(bases)
            nbp = length(basepairs)
            for n = 1:6
                @test sum(1 for _ in allseq(n; bases)) == numseq("."^n; nbases=nb)
            end
        end
    end

    @testset "allseq(dbn)" begin
        for (bases, basepairs) in TEST_BASE_BASEPAIRS
            for dbn in TEST_STRUCTURES_DBN
                nb = length(bases)
                nbp = length(basepairs)
                nseq = numseq(dbn; nbases=nb, nbasepairs=nbp)
                @test sum(1 for _ in allseq(dbn; bases, basepairs)) == nseq
                @test length(unique(collect(allseq(dbn; bases, basepairs)))) == nseq
            end
        end
    end

    @testset "allseq(pt)" begin
        for (bases, basepairs) in TEST_BASE_BASEPAIRS
            for dbn in TEST_STRUCTURES_DBN
                nb = length(bases)
                nbp = length(basepairs)
                pt = Pairtable(dbn)
                nseq = numseq(pt; nbases=nb, nbasepairs=nbp)
                @test sum(1 for _ in allseq(pt; bases, basepairs)) == nseq
                @test length(unique(collect(allseq(pt; bases, basepairs)))) == nseq
            end
        end
    end
end
