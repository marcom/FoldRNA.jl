
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

function design_random_ptarget(target::AbstractString, model;
                               niter::Integer=1000,
                               nbest::Integer=20)
    return design_random_ptarget(Pairtable(target), model; niter, nbest)
end
