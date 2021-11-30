using Test
using NucleicAcidFold: Alphabet, encode, decode

@testset "Alphabet" begin
    al = Alphabet("RNA", "ACGU")
    @test length(al) == 4
    @test encode(al, "GGCCUUAA") == [3, 3, 2, 2, 4, 4, 1, 1]
    @test decode(al, [3, 3, 2, 2, 4, 4, 1, 1]) == "GGCCUUAA"
    @test Alphabet("ABC") == Alphabet("", "ABC")
end
