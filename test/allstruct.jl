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
        end
    end
end
