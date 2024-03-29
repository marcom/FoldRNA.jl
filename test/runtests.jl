using Test

# show which testset is currently running
showtestset() = println(" "^(2 * Test.get_testset_depth()), "testing ",
                        Test.get_testset().description)

@testset verbose=true "FoldRNA" begin
    showtestset()
    include("fixedsize-priorityqueue.jl")

    include("alphabet.jl")
    include("semiring-log.jl")
    include("semiring-minplus.jl")

    include("pairtable.jl")
    include("allseq.jl")
    include("allstruct.jl")
    include("bpmodel.jl")
    include("loopmodel.jl")
    include("numstruct.jl")
    include("fold.jl")
    include("fold-bpmodel-score.jl")
    include("fold-loopmodel-score.jl")
    include("exhaustive.jl")
    include("design.jl")
    include("readparam-viennarna.jl")
end
