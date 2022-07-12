using Test
using FoldRNA: Alphabet, encode, decode

@testset "Alphabet" begin
    @test Alphabet("AB") == Alphabet("AB")
    al = Alphabet("RNA", "ACGU")
    @test length(al) == 4
    @test encode(al, "GGCCUUAA") == [3, 3, 2, 2, 4, 4, 1, 1]
    @test decode(al, [3, 3, 2, 2, 4, 4, 1, 1]) == "GGCCUUAA"
    @test Alphabet("ABC") == Alphabet("", "ABC")
    @test Alphabet(collect("ABC")) == Alphabet("", "ABC")
    @test encode(al, 'G') == [3]
    @test encode(al, ['A', 'G']) == [1, 3]
end
