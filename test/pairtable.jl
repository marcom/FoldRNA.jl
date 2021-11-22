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
            [ StrandInfo(1, 1, false), StrandInfo(2, 3, false), StrandInfo(4, 8, false) ]
        )
        @test length(pt) == 8
    end

end
