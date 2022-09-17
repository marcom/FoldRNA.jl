
function design_random_ptarget(target::Pairtable, model;
                               niter::Integer=1000, nbest::Integer=20)
    best = FixedsizePQ{String,Float64}(nbest)
    for i = 1:niter
        seq = randseq(target)
        seq in keys(best) && continue
        ptarget = prob_of_struct(Fold(seq, model), target)
        enqueue!(best, seq, ptarget)
    end
    return collect(best)
end

function design_random_ptarget(target::AbstractString, model; kwargs...)
    return design_random_ptarget(Pairtable(target), model; kwargs...)
end

function adaptivewalk(initial_x, loss, move!; maxiter=100, prob_accept_worse=0.01, show_trace::Bool=false)
    x = initial_x
    old_x = initial_x
    old_loss = loss(x)
    best_x = copy(old_x)
    best_loss = old_loss
    for iter = 1:maxiter
        new_x = copy(old_x)
        move!(new_x)
        new_loss = loss(new_x)
        if new_loss < old_loss
            best_loss = new_loss
            best_x = copy(new_x)
            if show_trace
                @show best_loss
            end
        end
        if new_loss < old_loss || rand() < prob_accept_worse
            old_x = new_x
            old_loss = new_loss
        end
    end
    return best_loss, best_x
end

function design_greedy_ptarget(target::Pairtable, model; startseq::AbstractString="", kwargs...)
    initial_x = collect(startseq == "" ? randseq(target) : startseq)
    loss(seq) = -prob_of_struct(Fold(join(seq), model), target)
    # TODO: hardcoded, filter alphabet from wildcards
    bases = ['A', 'C', 'G', 'U']
    # TODO: hardcoded, filter basepairs from wildcards
    #basepairs = [(bases[Tuple(I)[1]], bases[Tuple(I)[2]]) for i = 1:length(bases), j = 1:length(bases) if model.bptype[I] != 0]
    basepairs = [ ('A','U'), ('U','A'), ('C','G'), ('G','C'), ('G','U'), ('U','G')]
    function design_move!(seq::Vector{Char}, target::Pairtable)
        i = rand(firstindex(seq):lastindex(seq))
        if isbpopening(target, i) || isbpclosing(target, i)
            # change both pairs of a basepair
            j = target.pairs[i]
            ci, cj = rand(filter(bp -> bp != (seq[i],seq[j]), basepairs))
            seq[i] = ci
            seq[j] = cj
        else
            # change single base
            ci = rand(filter(c -> c != seq[i], bases))
            seq[i] = ci
        end
    end

    minus_p_target, aseq = adaptivewalk(initial_x, loss,
                                        seq -> design_move!(seq, target);
                                        prob_accept_worse=0.0,
                                        kwargs...)
    return -minus_p_target, join(aseq)
end

design_greedy_ptarget(target::AbstractString, model; kwargs...) =
    design_greedy_ptarget(Pairtable(target), model; kwargs...)
