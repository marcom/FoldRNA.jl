using Test
using NucleicAcidFold: Alphabet, encode, decode

@testset "Alphabet" begin
    al = Alphabet("RNA", "ACGU")
    @test encode(al, "GGCCUUAA") == [3, 3, 2, 2, 4, 4, 1, 1]
    @test decode(al, [3, 3, 2, 2, 4, 4, 1, 1]) == "GGCCUUAA"
end
