using Test
using FoldRNA: LogSR

@testset "Log semiring" begin
    showtestset()
    # TODO: test all properties from https://en.wikipedia.org/wiki/Semiring
    x = LogSR(1.0)
    @test x isa Number
    y = LogSR(2.0)
    @test (x + y) isa typeof(x)
    @test (x * y) isa typeof(x)
    @test (x / y) isa typeof(x)
    @test isapprox(float(y), exp(y.val))
    @test log(x) == x.val
    @test log(y) == y.val

    # test Base.show
    iobuf = IOBuffer()
    show(iobuf, MIME"text/plain"(), x)
    @test length(take!(iobuf)) > 0

    for T in (Float32, Float64, BigFloat)
        @test LogSR{T} <: Number
        s = LogSR{T}(T(4.2069))
        p = LogSR{T}(T(3.1415))
        z = zero(LogSR{T})
        o = one(LogSR{T})
        # test for type
        @test typeof(s) === typeof(p) == typeof(z) === typeof(o) === LogSR{T}
        for x in (s + p, p + s, s * p, p * s)
            @test x isa LogSR{T}
        end
        # + has to be commutative
        @test s + p == p + s
        # algebraic properties for one and zero
        for x in (s, p, z, o)
            @test x + z == x
            @test z + x == x
            @test x * z == z
            @test z * x == z
            @test x * o == x
            @test o * x == x
        end
        # +, * with ::T
        s = LogSR{T}(T(2.345))
        for x in (T(1.23), T(3.56))
            @test typeof(x + s) === LogSR{T}
            @test typeof(s + x) === LogSR{T}
            @test typeof(x * s) === LogSR{T}
            @test typeof(s * x) === LogSR{T}
        end
        # array broadcasting
        a = zeros(LogSR{T}, 10)
        a .= one(LogSR{T})
        @test a == ones(LogSR{T}, 10)
        a[:] .= zero(LogSR{T})
        @test a == zeros(LogSR{T}, 10)
    end
end
