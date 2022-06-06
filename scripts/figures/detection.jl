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

include(srcdir("main.jl"))

let conditioned_on = Rate(),
    postproc! = zscore!

n_signals = 1
boundary = Periodic()#Extended(50)
boundary_width = get_boundary_width(boundary)
n_trials=15
n_tests=15
α = 0.05 / 14

MEA_n = 100
MEA_t = 101
MEA_spikes = 267

raster_size = (MEA_n, (MEA_t*2)+(boundary_width*4))
signal_jitters = (12, 12)
lag_extents = (25, 25)
noise_rate = MEA_spikes / (MEA_n * MEA_t)

only_rate_raster = zeros(Bool, MEA_n, MEA_t + boundary_width*2)
@assert boundary_width == 0
only_rate_raster[1:MEA_spikes] .= true
precalced_postproc! = precalculate(postproc!, IndBernoulli(), conditioned_on, [only_rate_raster], boundary, lag_extents)

trial_results = @time run_peristimulus_tests(
    raster_size, boundary, lag_extents, 
    n_signals, signal_jitters, noise_rate,
    n_trials, n_tests,
    precalced_postproc!, α
)

plt = data(trial_results) * mapping(:signal_motif => "signal_motif", :detected_motif => "detected motif", color=:proportion_detected => "prop. tests positive") * visual(markersize=30)

axis = (; title = "$(n_trials) trials | $(n_tests) tests | $(n_signals) signals")
fig = draw(plt; axis)

outputdir = plotsdir("particular_trial_detects")
mkpath(outputdir)
save(joinpath(outputdir, "crossmotif_significances_$(fn2str(postproc!))|$(obj2str(conditioned_on))_trials$(n_trials)_tests$(n_tests)_signals$(n_signals)_size$(raster_size)_lags$(lag_extents).png"), fig)

end