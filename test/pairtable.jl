using Test
using NucleicAcidFold
using NucleicAcidFold: hasbp, isunpaired, isbpopening, isbpclosing

const TEST_DBN_SINGLESTRAND = [
    ".",
    ".....",
    "(.(..(..).(((...))))..)",
    "((..))..((...))",
    "...((..[[.))...((..]].))",
    ".(.[..{..).].}",
]
const TEST_DBN_MULTISTRAND = [
    ".+.",
    "(...).+(.)",
    ".[[.+.(..)...(((.]]..)))",
]
const TEST_DBN = vcat(TEST_DBN_SINGLESTRAND, TEST_DBN_MULTISTRAND)

@testset "Pairtable" begin
    @testset "constructor, ==, length" begin
        pt = Pairtable(4)
        @test pt == Pairtable(zeros(Int, 4), [StrandInfo(1, 4)])
        @test length(pt) == 4

        pt = Pairtable(2, 3)
        @test pt == Pairtable(zeros(Int, 5), [StrandInfo(1, 2),
                                              StrandInfo(3, 5)])
        @test length(pt) == 5

        pt = Pairtable(1, 2, 5)
        @test pt == Pairtable(zeros(Int, 8), [StrandInfo(1, 1),
                                              StrandInfo(2, 3),
                                              StrandInfo(4, 8)])
        @test length(pt) == 8
    end

    @testset "String(::Pairtable)" begin
        @test Pairtable([0, 3, 2, 0], [StrandInfo(1,4)]) |> String == ".()."
        @test Pairtable([0, 4, 0, 2, 0], [StrandInfo(1,2), StrandInfo(3,5)]) |> String == ".(+.)."
    end

    @testset "Parse non-pseudoknotted single" begin
        @test Pairtable(".")     == Pairtable([0], [StrandInfo(1, 1)])
        @test Pairtable("(.)")   == Pairtable([3, 0, 1], [StrandInfo(1, 3)])
        @test Pairtable("(.())") == Pairtable([5, 0, 4, 3, 1],
                                              [StrandInfo(1, 5)])
    end

    @testset "Parse non-pseudoknotted complex" begin
        @test Pairtable("(+)") == Pairtable([2, 1], [StrandInfo(1, 1),
                                                     StrandInfo(2, 2)])
        @test Pairtable("(().)+()") ==
            Pairtable([5, 3, 2, 0, 1, 7, 6], [StrandInfo(1, 5),
                                              StrandInfo(6, 7)])
    end

    @testset "Parse pseudoknotted single" begin
        @test Pairtable("([)]") == Pairtable([3, 4, 1, 2], [StrandInfo(1, 4)])
    end

    @testset "Parse pseudoknotted complex" begin
        @test Pairtable("(+[+)]") ==
            Pairtable([3, 4, 1, 2], [StrandInfo(1, 1), StrandInfo(2, 2),
                                     StrandInfo(3, 4)])
    end

    @testset "Parse illegal structures" begin
        @test_throws ErrorException Pairtable("")
        @test_throws ErrorException Pairtable(")")
        @test_throws ErrorException Pairtable("(")
        @test_throws ErrorException Pairtable("())")
        @test_throws ErrorException Pairtable("()(")
        @test_throws ErrorException Pairtable(")()")
        @test_throws ErrorException Pairtable("(()")
        @test_throws ErrorException Pairtable("[)")
        @test_throws ErrorException Pairtable("(])")
        @test_throws ErrorException Pairtable("()a")
        @test_throws ErrorException Pairtable("+")
        @test_throws ErrorException Pairtable("+(")
        @test_throws ErrorException Pairtable("+.")
        @test_throws ErrorException Pairtable(".+")
    end

    @testset "hash" begin
        @test typeof(hash(Pairtable("(...)"))) == UInt
        for s in TEST_DBN
            @test hash(Pairtable(s)) == hash(Pairtable(s))
        end
    end

    @testset "hasbp" begin
        @test hasbp(Pairtable(".(...)"), 2, 6)
        @test !hasbp(Pairtable("....."), 2, 5)
    end

    @testset "isunpaired" begin
        pt = Pairtable(".(...).")
        @test isunpaired(pt, 1)
        @test !isunpaired(pt, 2)
        @test isunpaired(pt, 3)
        @test isunpaired(pt, 4)
        @test isunpaired(pt, 5)
        @test !isunpaired(pt, 6)
        @test isunpaired(pt, 7)
    end

    @testset "isbpopening" begin
        pt = Pairtable(".((.))()")
        @test !isbpopening(pt, 1)
        @test isbpopening(pt, 2)
        @test isbpopening(pt, 3)
        @test !isbpopening(pt, 4)
        @test !isbpopening(pt, 5)
        @test !isbpopening(pt, 6)
        @test isbpopening(pt, 7)
        @test !isbpopening(pt, 8)
    end

    @testset "isbpclosing" begin
        pt = Pairtable("(.)().(())")
        @test !isbpclosing(pt, 1)
        @test !isbpclosing(pt, 2)
        @test isbpclosing(pt, 3)
        @test !isbpclosing(pt, 4)
        @test isbpclosing(pt, 5)
        @test !isbpclosing(pt, 6)
        @test !isbpclosing(pt, 7)
        @test !isbpclosing(pt, 8)
        @test isbpclosing(pt, 9)
        @test isbpclosing(pt, 10)
    end

    @testset "randseq" begin
        # length tests
        for n = 1:10
            @test length(randseq(n)) == n
            for bases in ["AB", "ACG"]
                @test length(randseq(n; bases)) == n
            end
        end
        for dbn in TEST_DBN_SINGLESTRAND
            pt = Pairtable(dbn)
            @test length(randseq(pt)) == length(randseq(dbn)) == length(pt)
            for (bases, basepairs) in [
                "AB" => ["AB"],
                "ACG" => ["CG", "GC"],
                "ACGUP" => ["AU", "UA", "AP", "PA", "CG", "GC", "GU", "UG", "GP", "PG"], ]
                @test length(randseq(pt; bases, basepairs)) == length(randseq(dbn; bases, basepairs)) == length(pt)
            end
        end
        # test random sequence generated
        bases = "ABC"
        basepairs = ["XY", "ZW",]
        for dbn in TEST_DBN_SINGLESTRAND
            pt = Pairtable(dbn)
            for _ = 1:5
                seq = randseq(pt; bases, basepairs)
                for i = 1:length(pt)
                    if isunpaired(pt, i)
                        @test seq[i] in bases
                    elseif isbpopening(pt, i)
                        @test seq[i] * seq[pt.pairs[i]] in basepairs
                    end
                end
            end
        end
    end
end
