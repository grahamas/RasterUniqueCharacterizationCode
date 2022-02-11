using DrWatson; quickactivate(@__DIR__)

using CSV, Tables, IterTools
using DataFrames
using LoopVectorization
using CairoMakie
using OffsetArrays
using Dates
using DataStructures
using TripleCorrelations

include(srcdir("rasters.jl"))
include(srcdir("motif_classes.jl"))
#include(srcdir("triple_correlation.jl"))
include(srcdir("latex_table.jl"))
include(srcdir("plotting.jl"))
include(srcdir("lag_motifs.jl"))
