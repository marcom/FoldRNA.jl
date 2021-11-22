using Test
using NucleicAcidFold

@testset "Pairtable" begin
    @testset "constructor, ==, length" begin
        pt = Pairtable(4)
        @test pt == Pairtable(zeros(Int, 4), [StrandInfo(1, 4, false)])
        @test length(pt) == 4

        pt = Pairtable(2, 3)
        @test pt == Pairtable(
            zeros(Int, 5),
            [StrandInfo(1, 2, false), StrandInfo(3, 5, false)]
        )
        @test length(pt) == 5

        pt = Pairtable(1, 2, 5)
        @test pt == Pairtable(
            zeros(Int, 8),
            [StrandInfo(1, 1, false),
             StrandInfo(2, 3, false),
             StrandInfo(4, 8, false)]
        )
        @test length(pt) == 8
    end

    @testset "Parse non-pseudoknotted single" begin
        @test Pairtable(".")     == Pairtable([0], [StrandInfo(1, 1, false)])
        @test Pairtable("(.)")   == Pairtable([3, 0, 1],
                                              [StrandInfo(1, 3, false)])
        @test Pairtable("(.())") == Pairtable([5, 0, 4, 3, 1],
                                              [StrandInfo(1, 5, false)])
    end

    @testset "Parse non-pseudoknotted complex" begin
        @test Pairtable("(+)") ==
            Pairtable([2, 1],
                      [StrandInfo(1, 1, false), StrandInfo(2, 2, false)])
        @test Pairtable("(().)+()") ==
            Pairtable([5, 3, 2, 0, 1, 7, 6],
                      [StrandInfo(1, 5, false), StrandInfo(6, 7, false)])
    end

    @testset "Parse pseudoknotted single" begin
        @test Pairtable("([)]") ==
            Pairtable([3, 4, 1, 2], [StrandInfo(1, 4, false)])
    end

    @testset "Parse pseudoknotted complex" begin
        @test Pairtable("(+[+)]") ==
            Pairtable([3, 4, 1, 2],
                      [StrandInfo(1, 1, false),
                       StrandInfo(2, 2, false),
                       StrandInfo(3, 4, false)])
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
end
