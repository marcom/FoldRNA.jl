using Test
using FoldRNA: MinPlusSR

@testset "MinPlus semiring" begin
    showtestset()
    # TODO: test all properties from https://en.wikipedia.org/wiki/Semiring
    x = MinPlusSR(1.0)
    @test x isa Number
    y = MinPlusSR(2.0)
    @test (x + y) isa typeof(x)
    @test (x * y) isa typeof(x)

    for T in (Float32, Float64, BigFloat)
        @test MinPlusSR{T} <: Number
        s = MinPlusSR{T}(T(4.2069))
        p = MinPlusSR{T}(T(3.1415))
        z = zero(MinPlusSR{T})
        o = one(MinPlusSR{T})
        # test for type
        @test typeof(s) === typeof(p) == typeof(z) === typeof(o) === MinPlusSR{T}
        for x in (s + p, p + s, s * p, p * s)
            @test x isa MinPlusSR{T}
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
        s = MinPlusSR{T}(T(2.345))
        for x in (T(1.23), T(3.56))
            @test typeof(x + s) === MinPlusSR{T}
            @test typeof(s + x) === MinPlusSR{T}
            @test typeof(x * s) === MinPlusSR{T}
            @test typeof(s * x) === MinPlusSR{T}
        end
        # array broadcasting
        a = zeros(MinPlusSR{T}, 10)
        a .= one(MinPlusSR{T})
        @test a == ones(MinPlusSR{T}, 10)
        a[:] .= zero(MinPlusSR{T})
        @test a == zeros(MinPlusSR{T}, 10)
    end
end
