using Test
using NucleicAcidFold: allstruct

@testset "allstruct" begin
    for hpmin = 0:3
        for n = 1:12
            @test length(unique(allstruct(n; hpmin))) == numstruct(n; hpmin)
            let c = 0
                allstruct(n; hpmin) do pt c += 1 end
                @test c == numstruct(n; hpmin)
            end
        end

        n = 12
        for _ = 1:5
            seq = randseq(n)
            @test length(unique(allstruct(seq; hpmin))) == numstruct(seq; hpmin)
            let c = 0
                allstruct(seq; hpmin) do pt c += 1 end
                @test c == numstruct(seq; hpmin)
            end

            seq = randseq(n; bases="ACGUN")
            @test length(unique(allstruct(seq; hpmin))) == numstruct(seq; hpmin)

            # test on Vector{Int} instead of String
            bases = "XYZW"  # deliberately different alphabet
            seq = randseq(n; bases)
            allbp = [(2,4), (4,2), (1,3), (3,1), (1,4), (4,1)]
            enc_seq = map(c -> findfirst(c, bases), collect(seq)) :: Vector{Int}
            @test length(unique(allstruct(enc_seq; hpmin, canbp=(s,i,j)->(s[i],s[j]) ∈ allbp))) ==
                numstruct(enc_seq; hpmin) do s,i,j return (s[i],s[j]) ∈ allbp end
        end
    end
end
