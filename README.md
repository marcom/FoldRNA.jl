# FoldRNA

This is an implementation of nucleic acid secondary structure
prediction and analysis algorithms in Julia.

## Installation

This package is not yet in the general Julia package registry. To
install it, press `]` in the Julia REPL to enter package mode and
install directly from the github repo with

```
add https://github.com/marcom/FoldRNA.jl
```

This will mean your installation directly tracks the `main` branch on
github.


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
- [Rfam.jl](https://github.com/cossio/Rfam.jl)
- [RNAstructure.jl](https://github.com/marcom/RNAstructure.jl)
- [ViennaRNA.jl](https://github.com/marcom/ViennaRNA.jl)
