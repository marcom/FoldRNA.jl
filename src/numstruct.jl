# References for number of secondary structures calculation:
#    [1] Waterman, Smith (1978) RNA Secondary Structure: A Complete
#        Mathematical Analysis. Math. Biosci., 42, 257–266.
#        (recursion formula, free energy minimization)
#    [2] Waterman (1978) Secondary structure of single-stranded
#        nucleic acids. In: Rota (ed.), Studies in Foundations and
#        Combinatorics number 1 in Advances in Mathematics, Academic
#        Press, pp. 167-212.
#        (proof of recursion)
#    [3] Stein, Waterman (1978) On some new sequences generalizing the
#        Catalan and Motzkin numbers. Discrete Math., 26, 261-272.
#        (contains table on p. 264 with values for n = 0:20, h = 0:6)

export numstruct
using OffsetArrays

"""
    numstruct(n::Integer; hpmin=3)

Number of unpseudoknotted secondary structures for a sequence of
length `n`.  The minimum allowed length of a hairpin loop is `hpmin`.
This function assumes that any base can base-pair with any other base.
"""
numstruct(n::Integer; hpmin::Integer=3) = numstruct(BigInt, n; hpmin)

function numstruct(::Type{T}, n::Integer; hpmin::Integer=3) where {T}
    n ≥ 0 || throw(ArgumentError("n must be non-negative (n = $n)"))
    hpmin ≥ 0 || throw(
        ArgumentError("hpmin must be non-negative (hpmin = $hpmin)"))
    n == 0 && return zero(T)
    # A[i]: number of structures for sequences of length i
    A = OffsetVector(ones(T, n+1), 0:n)
    for len = 2:n
        # case: first base is unpaired
        A[len] = A[len-1]
        # case: first base is paired, there are k bases contained
        # inside the base pair and len-k-2 bases outside the base pair
        for k = hpmin:len-2
            A[len] += A[k] * A[len - k - 2]
        end
    end
    return A[n]
end

"""
    numstruct(bp::Function, seq; hpmin=3)
    numstruct(T::Type, bp::Function, seq; hpmin=3)

Number of unpseudoknotted secondary structures for a sequence `seq`,
accumulated as type `T` (default is `BigInt`).  The function `bp(seq,
i, j)` indicates a multiplicative contribution to all the structures
that include the base-pair `(seq[i], seq[j])`.  The minimum allowed
length of a hairpin loop is `hpmin`."""

numstruct(bp::Function, seq; hpmin::Integer=3) =
    numstruct(BigInt, bp, seq; hpmin)

function numstruct(::Type{T}, bp::Function, seq; hpmin::Integer=3) where {T}
    length(seq) == 0 && return zero(T)
    return bpmodel(T, seq; hpmin, bp)[1, length(seq)]
end

"""
    numstruct(seq::AbstractString; hpmin=3)

Number of unpseudoknotted secondary structures for a sequence `seq`
given as an `AbstractString`.  The permissible base-pairs are
`$DEFAULT_BASEPAIRS`, and additionally the base `N` can base-pair with
any other base. The minimum allowed length of a hairpin loop is
`hpmin`. """
numstruct(seq::AbstractString; hpmin::Integer=3) =
    numstruct(default_canbp, seq; hpmin)
