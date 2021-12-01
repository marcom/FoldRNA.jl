# FixedsizePQ: fixed-size priority queue

# TODO: this should probably be put into DataStructures.jl

using DataStructures: PriorityQueue
import DataStructures: enqueue!, peek, dequeue!

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

peek(f::FixedsizePQ) = peek(f.pq)

function enqueue!(f::FixedsizePQ{K,V}, k::K, v::V) where {K,V}
    if length(f) < f.n
        enqueue!(f.pq, k, v)
    else
        ek, ev = peek(f.pq)
        if Base.Order.lt(f.pq.o, ev, v)
            dequeue!(f)
            enqueue!(f.pq, k, v)
        end
    end
    return f
end

enqueue!(f::FixedsizePQ{K,V}, p::Pair{K,V}) where {K,V} =
    enqueue!(f, first(p), last(p))

dequeue!(f::FixedsizePQ) = dequeue!(f.pq)
