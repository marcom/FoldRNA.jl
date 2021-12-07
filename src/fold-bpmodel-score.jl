# score, score_exp for Fold{BpModel}

struct Basepair
    i :: Int
    j :: Int
end

function score(f::Fold{M}, bp::Basepair) where {M <: BpModel}
    si = f.seq[bp.i]
    sj = f.seq[bp.j]
    BP = f.model.score_bp
    return BP[si, sj]
end

function score_exp(f::Fold{M}, bp::Basepair) where {M <: BpModel}
    si = f.seq[bp.i]
    sj = f.seq[bp.j]
    BP = f.model.score_bp
    RT = ustrip(f.model.RT)
    return exp(- BP[si, sj] / RT)
end
