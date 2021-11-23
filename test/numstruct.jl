using NucleicAcidFold
using Test, OffsetArrays

@testset "numstruct" begin
    @testset "numstruct(n::Int; hpmin)" begin
        # table[n,h] = number of unpseudoknotted secondary structures for
        #              sequence length n and minimum hairpin length h
        # rows start at n = 1
        # columns start at h = 0
        #
        # Table taken from
        #     Stein, Waterman (1978) On some new sequences generalizing the
        #     Catalan and Motzkin numbers. Discrete Math., 26, 261-272.
        #     (contains table on p. 264 with values for n = 0:20, h = 0:6)
        # Note: entry for n = 13, h = 0 is hard to read in online pdf, but
        # copy-and-pasting the text gives 41835.
        table = OffsetArray([
             #   h=0     h=1    h=2    h=3   h=4   h=5  h=6
                   1       1      1      1     1     1    1    # n =  1
                   2       1      1      1     1     1    1    # n =  2
                   4       2      1      1     1     1    1    # n =  3
                   9       4      2      1     1     1    1    # n =  4
                  21       8      4      2     1     1    1    # n =  5
                  51      17      8      4     2     1    1    # n =  6
                 127      37     16      8     4     2    1    # n =  7
                 323      82     33     16     8     4    2    # n =  8
                 835     185     69     32    16     8    4    # n =  9
                2188     423    146     65    32    16    8    # n = 10
                5798     978    312    133    64    32   16    # n = 11
               15511    2283    673    274   129    64   32    # n = 12
               41835    5373   1463    568   261   128   64    # n = 13
              113634   12735   3202   1184   530   257  128    # n = 14
              310572   30372   7050   2481  1080   517  256    # n = 15
              853467   72832  15605   5223  2208  1042  513    # n = 16
             2356779  175502  34705  11042  4528  2104 1029    # n = 17
             6536382  424748  77511  23434  9313  4256 2066    # n = 18
            18199284 1032004 173779  49908 19207  8624 4152    # n = 19
            50852019 2516347 390966 106633 39714 17504 8352    # n = 20
        ], 1:20, 0:6)
        for hpmin in axes(table,2)
            for n in axes(table,1)
                @test numstruct(n; hpmin) == table[n, hpmin]
                @test numstruct(BigInt, n; hpmin) == table[n, hpmin]
                @test numstruct(Int64, n; hpmin) == table[n, hpmin]
            end
        end
    end
end
