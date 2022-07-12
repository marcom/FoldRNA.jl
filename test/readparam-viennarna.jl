using Test
using FoldRNA: readparam_viennarna
import ViennaRNA_jll

@testset "readparam" begin
    @testset "viennarna" begin
        paramfile = joinpath(ViennaRNA_jll.artifact_dir, "share", "ViennaRNA", "rna_turner1999.par")
        dG, dH = readparam_viennarna(paramfile)
        @test dG isa LoopModel
        @test dH isa LoopModel
    end
end
