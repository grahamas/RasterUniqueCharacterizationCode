@quickactivate "RasterUniqueCharacterizationCode"

include(srcdir("load_results.jl"))
using CairoMakie; plot_ext = "png"
using AlgebraOfGraphics, Colors
using AlgebraOfGraphics: categorical


results_dir = joinpath(projectdir(), "remote_plots", "rate_divide_35trials_7sigs_MEA_100x101_261spk_periodic_2022_05_31-141003")
output_dir = joinpath(results_dir, "posthoc")
mkpath(output_dir)

results = load_results_df(joinpath(results_dir, "results.jld2"))

trials = 16
tests = results[1, :n_hypothesis_tests]
plt = data(filter(x -> x.n_trials == trials, results)) * mapping(:signal_motif => "signal motif", :detect_motif => "detected motif", color= :proportion_rejected => "prop. positive detections") * visual(markersize=30)

axis = (; title = "$(trials) trials | $(n_hypothesis_tests) tests | ")
fig = draw(plt; axis)