# score, score_exp for FoldPseq{BpModel}

function score(f::FoldPseq{T,M}, bp::Basepair) where {T, M <: BpModel}
    p = f.pseq
    B = f.model.score_bp
    return dot(p[:, bp.i], B, p[:, bp.j])
end

function score_exp(f::FoldPseq{T,M}, bp::Basepair) where {T, M <: BpModel}
    p = f.pseq
    B = f.model.score_bp
    RT = ustrip(f.model.RT)
    return exp(- dot(p[:, bp.i], B, p[:, bp.j]) / RT)
end
