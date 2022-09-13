# Extra functions that don't dispatch on Fold{LoopModel}, but only LoopModel

# TODO
# - put in loopmodel-param.jl
# - algorithms go in loopmodel-fold.jl
# - score_*(fold::Fold, ...) functions should call
#   the score_*(model::LoopModel, ...) functions

function score_hairpin_mismatch(model::LoopModel{T}, bptype_ij::Integer, bt_ip1::Integer,
                                bt_jm1::Integer, len::Integer) where {T}
    # TODO: hardcoded 3, make it part of LoopModel ?
    s = zero(T)
    if len == 3
        s += score_terminal_nonGC(model, bptype_ij)
    elseif len > 3
        s += score_mismatch(model, bptype_ij, bt_ip1, bt_jm1, model.mismatch_hairpin)
    end
    return s
end

function score_terminal_nonGC(model::LoopModel{T}, bptype::Integer) where {T}
    return model.terminal_nonGC_bp[bptype]
end

function score_mismatch(model::LoopModel{T}, bptype_ij::Integer, bt_k::Integer,
                        bt_l::Integer, mismatch_param) where {T}
    return mismatch_param[bptype_ij, bt_k, bt_l]
end
