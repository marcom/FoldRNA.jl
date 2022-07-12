using Test
using FoldRNA: design_random_ptarget

@testset "design" begin
    @testset "for BpModel" begin
        @testset "random_ptarget" begin
            model = RNA_BPMODEL
            target_dbn = "((...))"
            @test design_random_ptarget(target_dbn, model) isa Vector{Pair{String,Float64}}
        end
    end
end
