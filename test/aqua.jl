import Aqua
using FoldRNA

@testset "Aqua.test_all" begin
    showtestset()
    Aqua.test_all(FoldRNA)
end
