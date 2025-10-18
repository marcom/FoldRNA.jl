# FixedsizePQ: fixed-size priority queue

# TODO: this should probably be put into DataStructures.jl

using DataStructures: PriorityQueue

struct FixedsizePQ{K,V,Order}
    n :: Int
    pq :: PriorityQueue{K,V,Order}
end

function FixedsizePQ{K,V}(capacity::Integer, ord::Base.Order.Ordering=Base.Order.ForwardOrdering()) where {K,V}
    capacity <= 0 && throw(ArgumentError("capacity must be > 0"))
    return FixedsizePQ(capacity, PriorityQueue{K,V,typeof(ord)}(ord))
end

# TODO: Base.show

Base.length(f::FixedsizePQ) = length(f.pq)
Base.collect(f::FixedsizePQ) = collect(f.pq)
Base.keys(f::FixedsizePQ) = keys(f.pq)
Base.values(f::FixedsizePQ) = values(f.pq)
Base.first(f::FixedsizePQ) = first(f.pq)

function Base.push!(f::FixedsizePQ{K,V}, kv::Pair{K,V}) where {K,V}
    if length(f) < f.n
        push!(f.pq, kv)
    else
        # TODO: assert that length(f) == f.n
        ek, ev = first(f.pq)
        if Base.Order.lt(f.pq.o, ev, last(kv))
            popfirst!(f)
            push!(f.pq, kv)
        end
    end
    return f
end

Base.popfirst!(f::FixedsizePQ) = popfirst!(f.pq)
