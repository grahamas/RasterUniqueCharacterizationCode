# RasterUniqueCharacterizationCode

All figures and results in Deshpande, Smith, & van Drongelen (2022, submitted) were created using [Matlab (Mathworks)](https://www.mathworks.com/products/matlab.html) R2021b, and can be replicated by running the appropriately named script in the `matlab` directory.

### Julia

We also include Julia implementations of the tools we used for our analyses.

Table S2 was generated using the [Julia](https://julialang.org/) script `generate_all_lag_motifs_for_table.jl` (sans diagrams, which were created manually). 

Much of the Matlab code is replicated in Julia, though no Julia code was used to make any of the figures that appear in any published manuscript. To use the Julia part of this codebase (assuming [Julia v1.5+ is already installed](https://julialang.org/downloads/)):

0. Download this code base (e.g. via `git clone`).
1. Open a Julia console and enter these commands:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the Julia scripts and
everything should work out of the box, including correctly finding local paths.

### Prior versions

For the 2021 [preprint](https://www.biorxiv.org/content/10.1101/2021.08.16.456546v1), see the commit tagged [`preprint-v2`](https://github.com/grahamas/RasterUniqueCharacterizationCode/tree/preprint-v2). 

For completeness, we include the tag [`preprint-v1`](https://github.com/grahamas/RasterUniqueCharacterizationCode/tree/preprint-v1), which links to the code used in the first preprint submission (2021 Aug 17).
