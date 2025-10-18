using ResumableFunctions

"""
    allmoves(seq, pairtable)

Iterate over all possible one-basepair moves (insertions or deletions
of basepairs) for a sequence `seq` and structure `pairtable`.

Moves are yielded as tuples `(i, j)`.  If `i` and `j` are positive,
the move is a basepair insertion.  If `i` and `j` are negative, it is
a basepair deletion.
"""
@resumable function allmoves(seq::AbstractString, pt::Pairtable; hpmin::Int=3,
                             canbp::Function=default_canbp) :: Tuple{Int, Int}
    n = length(seq)
    if n != length(pt)
        throw(ArgumentError("Sequence and Pairtable must have same length."))
    end
    enclosing_bp = [(0, n+1)]
    sizehint!(enclosing_bp, (n÷2)+1)  # avoid extra allocations on push!() and pop!()
    for i = 1:n
        if isbpopening(pt, i)
            push!(enclosing_bp, (i, pt.pairs[i]))
            @yield (-i, -pt.pairs[i])
        elseif isbpclosing(pt, i)
            pop!(enclosing_bp)
        elseif isunpaired(pt, i)
            enclosing_j = enclosing_bp[end][2]
            k = i + 1
            while k < enclosing_j
                if isbpopening(pt, k)
                    k = pt.pairs[k] + 1
                elseif isunpaired(pt, k)
                    if k - i - 1 >= hpmin && canbp(seq, i, k)
                        @yield (i, k)
                    end
                    k += 1
                else
                    error("Illegal pairtable, position $k, enclosing_j = $enclosing_j, unexpected bp closing: $pt")
                end
            end
        else
            error("Illegal pairtable: $pt")
        end
    end
end

"""
    allneighbours(seq, pairtable)

For sequence `seq`, iterate over all neighbouring structures (one
basepair insertion or deletion away) from structure `pairtable`.
"""
@resumable function allneighbours(
        seq::AbstractString, pt::Pairtable; hpmin::Int=3,
        canbp::Function=default_canbp
    ) :: Pairtable
    for move in allmoves(seq, pt; hpmin, canbp)
        nb = deepcopy(pt)
        if move[1] < 0 && move[2] < 0
            # delete basepair (-i, -j)
            i, j = -move[1], -move[2]
            if i >= j
                error("Illegal move: $move")
            end
            if nb.pairs[i] != j || nb.pairs[j] != i
                error("Basepair $move not present in pairtable: $nb")
            end
            nb.pairs[i] = UNPAIRED
            nb.pairs[j] = UNPAIRED
        elseif move[1] > 0 && move[2] > 0
            # insert basepair (i, j)
            i, j = move
            if i >= j
                error("Illegal move: $move")
            end
            nb.pairs[i] = j
            nb.pairs[j] = i
        else
            error("Illegal move: $move")
        end
        @yield nb
    end
end

allneighbours(pt::Pairtable; hpmin::Int=3) = allneighbours("N"^length(pt), pt; hpmin)
