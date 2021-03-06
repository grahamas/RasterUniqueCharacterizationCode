@quickactivate "RasterUniqueCharacterizationCode"

using TripleCorrelations
using HypothesisTests, Random, Statistics
using ProgressMeter
using CairoMakie; plot_ext = "png"
using Dates
using Pkg
using Base.Threads

include(srcdir("roman_encode.jl"))
include(srcdir("peristimulus_testing.jl"))

# if !@isdefined(an_timeseries_dict) || force_redef
#     an_timeseries_dict = Dict()
# end

# if !@isdefined(peristimulus_an_results_dict) || force_redef
#     peristimulus_an_results_dict = Dict()
# end

force_redef = false

boundary = Periodic()
trials=1000
n_bootstraps=50
subdir = if boundary isa Periodic
    "AN_$(trials)trials_$(n_bootstraps)bs_periodic_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
elseif boundary isa ZeroPadded
    "AN_$(trials)trials_$(n_bootstraps)bs_zeropad_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
else
    error("Unrecognized boundary condition for TriCorr")
end
mkpath(plotsdir(subdir))
 
# @threads for motif_class_num = 2:14
# motif_class = roman_encode(motif_class_num)
#an_timeseries_dict[(motif_class,boundary)], peristimulus_an_results_dict[(motif_class,boundary)] = 
let n_size = 16, t_size = 60,
    n_max_jitter = 3, t_max_jitter = 2,
    n_lag = 6, t_lag = 5, t_step=2,
    t_window = 2t_lag + 1,
    noise_rate = 0.2;

# save_dir = plotsdir(subdir, "examples")
# mkpath(save_dir)
l_an_timeseries, trialavg_raster = an_timeseries_across_jittered_trials(motif_class_num, n_size, t_size, 1:n_size, -t_max_jitter:t_max_jitter, n_max_jitter, t_max_jitter, trials, noise_rate, boundary, (n_lag, t_lag), t_step, n_bootstraps; save_dir=false)
signal_raster = embedded_rand_motif(motif_class, n_size, t_size, -n_max_jitter:n_max_jitter, -t_max_jitter:t_max_jitter, n_max_jitter, t_max_jitter)

test_sizes = 1:max(trials÷10,1):trials
peristimulus_results = if haskey(peristimulus_an_results_dict, motif_class)
    peristimulus_an_results_dict[motif_class]
else
    l_motif_an_timeseries = [[contribs[motif_class_num] for contribs in timeseries] for timeseries ∈ l_an_timeseries]
    map(test_sizes) do test_size
        l_timeseries_sample = StatsBase.sample(l_motif_an_timeseries, test_size, replace=false)
        peristimulus_start, peristimulus_stop = calculate_jitter_peristimulus_window(t_max_jitter, t_window, t_step, t_size ÷ 2)
        test_peristimulus_difference(l_timeseries_sample, peristimulus_start, peristimulus_stop)
    end
end

noise_raster = rand(size(signal_raster)...) .< noise_rate
raster = min.(signal_raster .+ noise_raster, 1)

f_signal = heatmap(signal_raster', axis=(xlabel="time", ylabel="neuron"))
f_noise = heatmap(noise_raster', axis=(xlabel="time", ylabel="neuron"))
f_raster = heatmap(raster', axis=(xlabel="time", ylabel="neuron"))
f_trialavg_raster = heatmap(trialavg_raster', axis=(xlabel="time", ylabel="neuron"))

@show length(l_an_timeseries)
@show size(l_an_timeseries |> first)
f_motif_course = plot(1:t_step:(size(raster,2)-2t_lag), [a[motif_class_num] for a ∈ mean(l_an_timeseries)], axis=(xlabel="time", ylabel="avg motif $(motif_class) contrib"))
f_motif_control = plot(1:t_step:(size(raster,2)-2t_lag), [a[14] for a ∈ mean(l_an_timeseries)], axis=(xlabel="time", ylabel="avg motif XIV contrib"))

save(plotsdir(subdir,"signal_motif_$(motif_class)_AN.$(plot_ext)"), f_signal)
save(plotsdir(subdir,"noise_motif_$(motif_class)_AN.$(plot_ext)"), f_noise)
save(plotsdir(subdir,"raster_motif_$(motif_class)_AN.$(plot_ext)"), f_raster)
save(plotsdir(subdir,"trialavg_raster_motif_$(motif_class)_AN.$(plot_ext)"), f_trialavg_raster)

save(plotsdir(subdir,"motif_$(motif_class)_AN_timeseries_$(typeof(boundary)).$(plot_ext)"), f_motif_course)
save(plotsdir(subdir,"motif_XIV_given_$(motif_class)_AN_timeseries_$(typeof(boundary)).$(plot_ext)"), f_motif_control)

(l_an_timeseries, peristimulus_results)
end

# let l_an_timeseries=an_timeseries_dict[motif_class];
# end

