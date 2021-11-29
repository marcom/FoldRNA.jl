export energy, partfn
using OffsetArrays: OffsetMatrix
using Unitful: @u_str, Quantity, ustrip
using DataStructures: DefaultDict

struct BpModel{T}
    name :: String
    al :: Alphabet
    score_bp :: Matrix{T}
    RT :: Quantity
    unit :: Quantity
    hpmin :: Int
end
encode(m::BpModel, iter) = encode(m.al, iter)
decode(m::BpModel, iter) = decode(m.al, iter)

const RNA_BPMODEL = BpModel{Float64}(
    "RNA Nussinov-Jacobson model (basepairs)",
    Alphabet("RNA", "ACGU"),
    [ Inf Inf Inf -2.0; Inf Inf -3.0 Inf; Inf -3.0 Inf -1.0; -2.0 Inf -1.0 Inf ],
    RT37,
    1.0u"kcal/mol",
    DEFAULT_HPMIN
)

# Base-pair model recursions. O(n^3) time, O(n^2) space.
function bpmodel(::Type{T}, seq; hpmin::Integer=3, bp::Function) where {T}
    # TODO: check that indices of seq are 1:n
    n = length(seq)
    A = OffsetMatrix(Matrix{T}(undef, n+1, n+1), 0:n, 0:n)
    return bpmodel(A, seq; hpmin, bp)
end

function bpmodel(A::OffsetMatrix{T}, seq; hpmin::Integer=3, bp::Function) where {T}
    # TODO: check that indices of seq are 1:n
    n = length(seq)
    hpmin ≥ 0 || throw(ArgumentError("hpmin must be non-negative (hpmin = $hpmin)"))
    n > 0 || throw(ArgumentError("length(seq) must be ≥ 0 (length(seq) = $n)"))
    axes(A) == (0:n, 0:n) ||
        throw(ArgumentError("axes of A must be (0:n, 0:n), n = length(seq)"))
    A .= one(T)
    for d = 1:n
        for i = 1:n-d
            j = i + d
            # case: j is unpaired
            A[i,j] = A[i, j-1]
            # case: j is paired to k in [i, j-h-1]
            for k = i:j-hpmin-1
                A[i,j] += A[i, k-1] * A[k+1, j-1] * bp(seq, k, j)
            end
        end
    end
    return A
end
