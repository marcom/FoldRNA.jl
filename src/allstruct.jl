export allstruct
using ResumableFunctions

"""
    allstruct(seq::AbstractString; hpmin=3, canbp::Function=(s,i,j)->...)
    allstruct(n::Integer; hpmin=3, canbp::Function)
    allstruct(f::Function, seq::String; hpmin=3, canbp::Function)
    allstruct(f::Function, n::Integer; hpmin=3, canbp::Function)

Iterate over all unpseudoknotted secondary structures, either for a
given sequence `seq` or all structures of length `n`.  The methods
that take a function `f` as the first argument apply this function to
every structure.

Optional arguments: minimum hairpin size is `hpmin`, and the function
`canbp` specifies which bases can form basepairs.
"""
function allstruct(seq::AbstractString;
                   hpmin::Integer=3,
                   canbp::Function=default_canbp)
    return _allstruct_rf(seq; hpmin, canbp)
end

function allstruct(f::Function,
                   seq::AbstractString;
                   hpmin::Integer=3,
                   canbp::Function=default_canbp)
    return _allstruct_rf(f, seq; hpmin, canbp)
end

allstruct(n::Integer; hpmin::Integer=3) =
    allstruct("N"^n; hpmin, canbp=(s,i,j)->true)

allstruct(f::Function, n::Integer; hpmin::Integer=3) =
    allstruct(f, "N"^n; hpmin, canbp=(s,i,j)->true)

@resumable function _allstruct_rf(seq::AbstractString;
                                  hpmin::Integer, canbp::Function) :: Pairtable
    n = length(seq)
    pt = Pairtable(n)
    for s in _allstruct_rf_rec!(pt, seq, 1, n, hpmin, canbp)
        @yield deepcopy(pt)
    end
end

function _allstruct_rf(f, seq::AbstractString;
                       hpmin::Integer, canbp::Function)
    n = length(seq)
    pt = Pairtable(n)
    for s in _allstruct_rf_rec!(pt, seq, 1, n, hpmin, canbp)
        f(pt)
    end
end

@resumable function _allstruct_rf_rec!(pt::Pairtable, seq::AbstractString,
                                       i::Integer, n::Integer,
                                       hpmin::Integer, canbp::Function) :: Pairtable
    if i > n
        # empty base case
        @yield pt
        return
    end
    # case: i unpaired
    pt.pairs[i] = UNPAIRED
    for s in _allstruct_rf_rec!(pt, seq, i+1, n, hpmin, canbp)
        @yield pt
    end
    # case: i paired to k
    for k = i+hpmin+1:n
        if canbp(seq, i, k)
            pt.pairs[i] = k
            pt.pairs[k] = i
            for s1 in _allstruct_rf_rec!(pt, seq, i+1, k-1, hpmin, canbp)
                for s2 in _allstruct_rf_rec!(pt, seq, k+1, n, hpmin, canbp)
                    @yield pt
                end
            end
        end
    end
end
