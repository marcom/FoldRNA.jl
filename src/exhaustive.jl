function exhaustive_partfn(::Type{T}, seq::AbstractString, param::BPmodelParam;
                           hpmin::Integer=3) where {T}
    # TODO: take temperature from param
    Q = zero(T)
    for pt in allstruct(seq; hpmin)
        en = energy(seq, pt, param)
        Q += exp(- en / param.RT)
    end
    return Q
end

function exhaustive_partfn(seq::AbstractString, param::BPmodelParam;
                           hpmin::Integer=3)
    s = exhaustive_partfn(LogSR{Float64}, seq, param; hpmin)
    logQ = s.val
    return - param.RT * logQ
end
