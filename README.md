# RasterUniqueCharacterizationCode

## Matlab

This code base is using Matlab (Mathworks). To reproduce this project, do the following: 
Download triple_correlation_script

Figure 1 (not reproduced by this code)
Figure 2A: run code as is
Figure 2B: 
Figure 2C-D: Change noise scaling factor in Line 27 to 1 and signal scaling factor to 0, and run the code
Figure 2E-F: Comment Lines 30-39, uncomment lines 41-46, and run the code
Figure 2G-H: 

Figure 3A: Change frequency (f) in line 11 to f = 0.12, and run the code
Figure 3B: Uncomment line 19, comment line 20 for the phase shift, uncomment line 99, and run the code
Figure 3C: Change the noise scaling factor to 1, uncomment line 28, and run the code
Figure 3D: Change the signal scaling factor to 0, and run the code


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
