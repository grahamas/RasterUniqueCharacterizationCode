# RasterUniqueCharacterizationCode

All figures and results in Deshpande, Smith, & van Drongelen (2021, submitted; [preprint](https://www.biorxiv.org/content/10.1101/2021.08.16.456546v1)[^1]) were created using [Matlab (Mathworks)](https://www.mathworks.com/products/matlab.html) R2020b, and can be replicated by running `triple_correlation_script.mat` in the `scripts` directory.

   - Figure 1: no data (not reproduced by this code)
   - Figures 2-3: Run `triple_correlation_script.mat` in Matlab

### Julia

Table S2 was generated using the [Julia](https://julialang.org/) script `generate_all_lag_motifs_for_table.jl` (sans diagrams, which were created manually). 

Much of the Matlab code is replicated in Julia, though no Julia code was used to make any of the figures that appear in any published manuscript. To use the Julia part of this codebase (assuming [Julia v1.5+ is already installed](https://julialang.org/downloads/)):

0. Download this code base (e.g. via `git clone`).
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the Julia scripts and
everything should work out of the box, including correctly finding local paths.

   [^1]: Note that preprint v1 has figures differing from the included code, though displaying the same underlying data.
