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

function embedded_static_rand_motif(motif_class, n_size, t_size, n_max_jitter, t_max_jitter)
    raster = zeros(Bool, n_size, t_size)
    motif_coords = TripleCorrelations.rand_motif(motif_class, 0:0, 0:0, n_max_jitter, t_max_jitter)
    for motif_coord in motif_coords
        raster[CartesianIndex(motif_coord .+ (n_size ÷ 2, t_size ÷ 2))] = 1
    end
    return raster
end

if !@isdefined(detect_an_across_static_rand_trials) || force_redef
@memoize function detect_an_across_rand_static_trials(motif_class, n_size, t_size, n_max_jitter, t_max_jitter, trials, noise_rate, boundary, n_lag, t_lag, t_step, n_bootstraps)
    signal_raster = embedded_static_rand_motif(motif_class, n_size, t_size, n_max_jitter, t_max_jitter)
    l_an_timeseries, trialavg_raster = detect_an_across_trials(signal_raster, trials, noise_rate, boundary, n_lag, t_lag, t_step, n_bootstraps)
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

for motif_class_num = 2:5
motif_class = roman_encode(motif_class_num)
an_timeseries_dict[(motif_class,boundary)], peristimulus_an_results_dict[(motif_class,boundary)] = let n_size = 16, t_size = 60,
    n_max_jitter = 4, t_max_jitter = 2,
    n_lag = 6, t_lag = 5, t_step=2,
    t_window = 2t_lag + 1,
    noise_rate = 0.2;

l_an_timeseries, signal_raster, trialavg_raster = detect_an_across_rand_static_trials(motif_class, n_size, t_size, n_max_jitter, t_max_jitter, trials, noise_rate, boundary, n_lag, t_lag, t_step, n_bootstraps)

test_sizes = 1:max(trials÷10,1):trials
peristimulus_results = if haskey(peristimulus_an_results_dict, motif_class)
    peristimulus_an_results_dict[motif_class]
else
    l_motif_an_timeseries = [[contribs[motif_class_num] for contribs in timeseries] for timeseries ∈ l_an_timeseries]
    map(test_sizes) do test_size
        l_timeseries_sample = rand(l_motif_an_timeseries, test_size)
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

