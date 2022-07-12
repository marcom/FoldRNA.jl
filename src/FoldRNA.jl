module FoldRNA

include("fixedsize-priorityqueue.jl")

include("alphabet.jl")
include("semiring-log.jl")
include("semiring-minplus.jl")

include("defaults.jl")
include("pairtable.jl")
include("allseq.jl")
include("allstruct.jl")
include("bpmodel.jl")
include("loopmodel.jl")
include("numstruct.jl")
include("fold.jl")
include("fold-bpmodel-score.jl")
include("fold-loopmodel-score.jl")
include("exhaustive.jl")
include("design.jl")
include("readparam-viennarna.jl")

# LoopModel parameters from ViennaRNA
import ViennaRNA_jll
# only export free energy params for now, not the enthalpies
export RNA_TURNER1999, RNA_TURNER2004, RNA_ANDRONESCU2007, RNA_LANGDON2018

const RNA_TURNER1999, RNA_TURNER1999_DH = readparam_viennarna(
    joinpath(ViennaRNA_jll.artifact_dir, "share", "ViennaRNA", "rna_turner1999.par")
)
const RNA_TURNER2004, RNA_TURNER2004_DH = readparam_viennarna(
    joinpath(ViennaRNA_jll.artifact_dir, "share", "ViennaRNA", "rna_turner2004.par")
)
const RNA_ANDRONESCU2007, RNA_ANDRONESCU2007_DH = readparam_viennarna(
    joinpath(ViennaRNA_jll.artifact_dir, "share", "ViennaRNA", "rna_andronescu2007.par")
)
const RNA_LANGDON2018, RNA_LANGDON2018_DH = readparam_viennarna(
    joinpath(ViennaRNA_jll.artifact_dir, "share", "ViennaRNA", "rna_langdon2018.par")
)
# TODO: DNA parameters maybe have some issues? (see comments in
#       ViennaRNA dna_*.par files).  Also: are these DNA params
#       measured at 37 deg Celsius?

end # module FoldRNA
