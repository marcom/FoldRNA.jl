using Test
using FoldRNA: FixedsizePQ

@testset "FixedsizePQ" begin
    showtestset()
    fpq = FixedsizePQ{String,Float64}(2)
    @test length(fpq) == 0
    @test collect(fpq) == []
    @test length(keys(fpq)) == 0
    @test length(values(fpq)) == 0
    push!(fpq, "a" => 3.0)
    @test length(fpq) == 1
    @test first(fpq) == ("a" => 3.0)
    @test popfirst!(fpq) == ("a" => 3.0)
    @test length(fpq) == 0

    push!(fpq, "a" => 1.0)
    push!(fpq, "b" => 2.0)
    push!(fpq, "c" => 3.0)
    @test length(fpq) == 2
    @test Set(collect(fpq)) == Set(["b" => 2.0, "c" => 3.0])
    @test Set(keys(fpq)) == Set(["b", "c"])
    @test Set(values(fpq)) == Set([2.0, 3.0])
    @test first(fpq) == ("b" => 2.0)
    @test popfirst!(fpq) == ("b" => 2.0)
    @test length(fpq) == 1
    @test Set(collect(fpq)) == Set(["c" => 3.0])
    @test first(fpq) == ("c" => 3.0)
end
