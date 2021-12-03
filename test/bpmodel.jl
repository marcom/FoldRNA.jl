using Test
using NucleicAcidFold: bpmodel, bpmodel_bpp
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
    @testset "bpmodel_bpp" begin
        fold = Fold("GGGAAACCC", RNA_BPMODEL)
        n = length(fold)
        A = bpmodel(LogSR{Float64}, fold; hpmin=fold.model.hpmin,
                    bp = (f,i,j) -> 1.0)
        p = bpmodel_bpp(A, fold; hpmin=fold.model.hpmin,
                        bp = (f,i,j) -> 1.0)
        @test p isa Matrix
        @test size(p) == (n,n)
        @test all(x -> 0.0 <= x <= 1.0, p)
    end
end
