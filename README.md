# FoldRNA

[![Build Status](https://github.com/marcom/FoldRNA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/marcom/FoldRNA.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is an implementation of nucleic acid secondary structure
prediction and analysis algorithms in Julia.

## Installation

Press `]` in the Julia REPL to enter package mode and enter

```
add FoldRNA
```


## Basic usage

```julia
using FoldRNA

# Nussinov-Jacobson model (basepair model)

f = Fold("GGGAAACCC", RNA_BPMODEL)
mfe(f)
partfn(f)
bpp(f)
energy(f, "(((...)))")
prob_of_struct(f, "(((...)))")

# Nearest-neighbour energy model (loop model)
#
# Note: there are still some errors in the implementation of the
#       energy calculation in the loop decomposition, so mfe() and
#       partfn() give wrong results, but energy() should work

f = Fold("GGGAAACCC")    # defaults to RNA_TURNER2004 energy params
mfe(f)
partfn(f)
# bpp(f) not implemented yet for loopmodel
prob_of_struct(f, "(((...)))")
```

## Related Julia packages for RNA secondary structure

- [Infernal.jl](https://github.com/cossio/Infernal.jl)
- [LinearFold.jl](https://github.com/marcom/LinearFold.jl)
- [PlotRNA.jl](https://github.com/marcom/PlotRNA.jl)
- [Rfam.jl](https://github.com/cossio/Rfam.jl)
- [RNAstructure.jl](https://github.com/marcom/RNAstructure.jl)
- [ViennaRNA.jl](https://github.com/marcom/ViennaRNA.jl)
