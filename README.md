# RasterUniqueCharacterizationCode

## Matlab

All figures are made by running `triple_correlation_script.mat`.

   - Figure 1: no data (not reproduced by this code)
   - Figures 2-3: Run code as is

## Julia

Table S2 was generated using `generate_all_lag_motifs_for_table.jl` (sans diagrams, which were created manually). 

Much of the Matlab code is replicated in Julia, though not used to make any figures.

### Installing Julia codebase

This code base uses the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> RasterUniqueCharacterizationCode

To (locally) instantiate this julia code with its dependencies, do the following:

0. Download this code base (e.g. via `git clone`).
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the julia scripts and
everything should work out of the box, including correctly finding local paths.
