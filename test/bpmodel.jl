using Test
using NucleicAcidFold: bpmodel
using Unitful: @u_str

@testset "bpmodel" begin
    @testset "bpmodel" begin
        # check that number of visited structures is correct
        function mynumstruct(seq; hpmin)
            A = bpmodel(BigInt, seq;
                        hpmin, bp = (s,i,j) -> true)
            return A[1, length(seq)]
        end
        for n = 30:33
            for hpmin = 0:3
                @test mynumstruct("N"^n; hpmin) == numstruct(n; hpmin)
            end
        end
    end
end
