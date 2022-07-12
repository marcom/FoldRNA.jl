using Test
using FoldRNA: FixedsizePQ, enqueue!, dequeue!, peek

@testset "FixedsizePQ" begin
    fpq = FixedsizePQ{String,Float64}(2)
    @test length(fpq) == 0
    @test collect(fpq) == []
    @test length(keys(fpq)) == 0
    @test length(values(fpq)) == 0
    enqueue!(fpq, "a" => 3.0)
    @test length(fpq) == 1
    @test peek(fpq) == ("a" => 3.0)
    @test dequeue!(fpq) == "a"
    @test length(fpq) == 0

    enqueue!(fpq, "a" => 1.0)
    enqueue!(fpq, "b" => 2.0)
    enqueue!(fpq, "c" => 3.0)
    @test length(fpq) == 2
    @test Set(collect(fpq)) == Set(["b" => 2.0, "c" => 3.0])
    @test Set(keys(fpq)) == Set(["b", "c"])
    @test Set(values(fpq)) == Set([2.0, 3.0])
    @test peek(fpq) == ("b" => 2.0)
    @test dequeue!(fpq) == "b"
    @test length(fpq) == 1
    @test Set(collect(fpq)) == Set(["c" => 3.0])
    @test peek(fpq) == ("c" => 3.0)
end
