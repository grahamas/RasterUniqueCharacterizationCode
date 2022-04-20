@quickactivate "RasterUniqueCharacterizationCode"

using TripleCorrelations
using HypothesisTests, Random, Statistics
using ProgressMeter
using CairoMakie; plot_ext = "png"
using Dates
using Pkg

include(srcdir("roman_encode.jl"))
include(srcdir("peristimulus_testing.jl"))

if !@isdefined(an_timeseries_dict) || force_redef
    an_timeseries_dict = Dict()
end

if !@isdefined(peristimulus_an_results_dict) || force_redef
    peristimulus_an_results_dict = Dict()
end

if !@isdefined(detect_mid_an_across_trials) || force_redef
@memoize function detect_mid_an_across_trials(motif_class_num, n_pad, t_pad, n_reps, t_reps, trials, noise_rate, boundary, n_lag, t_lag, t_step, n_bootstraps; save_dir=false)
    motif_class = roman_encode(motif_class_num)
    signal_raster = TripleCorrelations.repeat_padded_motif(motif_class, n_pad, t_pad, n_reps, t_reps)
    l_an_timeseries, trialavg_raster = detect_an_across_trials(motif_class_num, signal_raster, trials,noise_rate, boundary, n_lag, t_lag, t_step, n_bootstraps; save_dir=save_dir)
    return (l_an_timeseries, signal_raster, trialavg_raster)
end
end
force_redef = false

boundary = Periodic()
trials=100
n_bootstraps=10
subdir = if boundary isa Periodic
    "AN_$(trials)trials_$(n_bootstraps)bs_periodic_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
elseif boundary isa ZeroPadded
    "AN_$(trials)trials_$(n_bootstraps)bs_zeropad_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
else
    error("Unrecognized boundary condition for TriCorr")
end
mkpath(plotsdir(subdir))

for motif_class_num = [5]
motif_class = roman_encode(motif_class_num)
an_timeseries_dict[(motif_class,boundary)], peristimulus_an_results_dict[(motif_class,boundary)] = let n_pad = 10, t_pad = 30, 
    n_reps = 1, t_reps = 1,
    n_lag = 6, t_lag = 5,
    trials=100, t_step=2,
    t_window = 2t_lag + 1,
    noise_rate = 0.2,
    n_bootstraps = 10;

save_dir = plotsdir(subdir, "examples")
mkpath(save_dir)
l_an_timeseries, signal_raster, trialavg_raster = detect_mid_an_across_trials(motif_class_num, n_pad, t_pad, n_reps, t_reps, trials, noise_rate, boundary, n_lag, t_lag, t_step, n_bootstraps; save_dir=save_dir)

test_sizes = 1:max(trials÷10,1):trials
peristimulus_results = if haskey(peristimulus_an_results_dict, motif_class)
    peristimulus_an_results_dict[motif_class]
else
    l_motif_an_timeseries = [[contribs[motif_class_num] for contribs in timeseries] for timeseries ∈ l_an_timeseries]
    map(test_sizes) do test_size
        l_timeseries_sample = StatsBase.sample(l_motif_an_timeseries, test_size, replace=false)
        motif_t_bounds = get_motif_t_bounds(motif_class)
        peristimulus_start, peristimulus_stop = calculate_peristimulus_window(t_pad, t_window, t_step, motif_t_bounds...)
        test_peristimulus_difference(l_timeseries_sample, peristimulus_start, peristimulus_stop)
    end
end

noise_raster = rand(size(signal_raster)...) .< noise_rate
raster = min.(signal_raster .+ noise_raster, 1)

f_signal = heatmap(signal_raster', axis=(xlabel="time", ylabel="neuron"))
f_noise = heatmap(noise_raster', axis=(xlabel="time", ylabel="neuron"))
f_raster = heatmap(raster', axis=(xlabel="time", ylabel="neuron"))
f_trialavg_raster = heatmap(trialavg_raster', axis=(xlabel="time", ylabel="neuron"))

f_motif_course = plot([a[motif_class_num] for a ∈ mean(l_an_timeseries)], axis=(xlabel="time", ylabel="avg motif $(motif_class) contrib"))
f_motif_control = plot([a[14] for a ∈ mean(l_an_timeseries)], axis=(xlabel="time", ylabel="avg motif XIV contrib"))

save(plotsdir(subdir,"signal_motif_$(motif_class)_AN.$(plot_ext)"), f_signal)
save(plotsdir(subdir,"noise_motif_$(motif_class)_AN.$(plot_ext)"), f_noise)
save(plotsdir(subdir,"raster_motif_$(motif_class)_AN.$(plot_ext)"), f_raster)
save(plotsdir(subdir,"trialavg_raster_motif_$(motif_class)_AN.$(plot_ext)"), f_trialavg_raster)

save(plotsdir(subdir,"motif_$(motif_class)_AN_timeseries_$(typeof(boundary)).$(plot_ext)"), f_motif_course)
save(plotsdir(subdir,"motif_XIV_given_$(motif_class)_AN_timeseries_$(typeof(boundary)).$(plot_ext)"), f_motif_control)

(l_an_timeseries, peristimulus_results)
end
end

# let l_an_timeseries=an_timeseries_dict[motif_class];
# end

