@quickactivate "RasterUniqueCharacterizationCode"

using TripleCorrelations
using HypothesisTests, Random, Statistics
using ProgressMeter
using CairoMakie; plot_ext = "png"
using Dates
using Pkg
using Base.Threads
using DataFrames, AlgebraOfGraphics, Colors
using JLD2

include(srcdir("roman_encode.jl"))
include(srcdir("peristimulus_testing.jl"))


raster_size = (150, 200)
lag_extents = (5,7)