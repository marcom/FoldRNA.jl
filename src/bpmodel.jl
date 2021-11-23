using OffsetArrays: OffsetMatrix

# Base-pair model recursions. O(n^3) time, O(n^2) space.
function bpmodel(T::Type, seq; hpmin::Integer=3, bp::Function)
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
