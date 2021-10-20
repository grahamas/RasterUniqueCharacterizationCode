# RasterUniqueCharacterizationCode

## Matlab

This code base is using Matlab (Mathworks). To reproduce this project, do the following: 

All figures are made by running `triple_correlation_script.mat`.

   - Figure 1: no data (not reproduced by this code)
   - Figures 2-3: Run code as is

## Julia

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> RasterUniqueCharacterizationCode

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
