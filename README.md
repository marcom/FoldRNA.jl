# FoldRNA

This is an experimental implementation of nucleic acid secondary
structure prediction and analysis algorithms in Julia.  Everything is
still in a state of flux and not guaranteed to work.

## Installation

This package is not yet in the general Julia package registry. To
install it, press `]` in the Julia REPL to enter package mode and
install directly from the github repo with

```
add https://github.com/marcom/FoldRNA.jl
```

This will mean your installation directly tracks the `main` branch on
github, and changes will be pulled directly from github every time you
update packages.  For more details, check out the
[section on "Managing Packages"](https://pkgdocs.julialang.org/v1/managing-packages/)
section of the Julia package manager docs.

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

## Related packages

If you are looking for more mature libraries for RNA folding in Julia,
please take a look at:

- [ViennaRNA.jl](https://github.com/marcom/ViennaRNA.jl)
- [LinearFold.jl](https://github.com/marcom/LinearFold.jl)
