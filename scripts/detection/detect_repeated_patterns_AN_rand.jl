@quickactivate "RasterUniqueCharacterizationCode"

using TripleCorrelations
using HypothesisTests, Random, Statistics
using ProgressMeter
using CairoMakie; plot_ext = "png"
using Dates
using Pkg
using Base.Threads
using DataFrames, AlgebraOfGraphics
using JLD2

include(srcdir("roman_encode.jl"))
include(srcdir("peristimulus_testing.jl"))

if !@isdefined(prior_results_dict) || force_redef
    prior_results_dict = Dict()
end
force_redef = false

boundary = Periodic()
trials=300
n_bootstraps=20
n_resamples=100
α = 0.05
results_key = (; boundary=boundary, trials=trials, n_bootstraps=n_bootstraps, n_resamples=n_resamples, α=α)
subdir = if boundary isa Periodic
    "AN_$(trials)trials_$(n_bootstraps)bs_periodic_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
elseif boundary isa ZeroPadded
    "AN_$(trials)trials_$(n_bootstraps)bs_zeropad_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
else
    error("Unrecognized boundary condition for TriCorr")
end
mkpath(plotsdir(subdir))

@threads for motif_class_num = 2:14
#for motif_class_num = 6:6
motif_class = offset_motif_numeral(motif_class_num)
prior_results_dict[merge((motif_class=motif_class,), results_key)] = 
let n_size = 16, t_size = 60,
    n_max_jitter = 3, t_max_jitter = 2,
    n_lag = 6, t_lag = 5, t_step=2,
    t_window = 2t_lag + 1,
    noise_rate = 0.2;

# save_dir = plotsdir(subdir, "examples")
# mkpath(save_dir)

test_sizes = 1:max(trials÷10,1):trials
l_an_timeseries, trialavg_raster, peristimulus_results = getkey(prior_results_dict, merge((motif_class=motif_class,), results_key), (begin
    l_an_timeseries, trialavg_raster = detect_an_across_jittered_trials(
        motif_class_num, n_size, t_size, 1:n_size, -t_max_jitter:t_max_jitter, n_max_jitter, t_max_jitter, trials, noise_rate, boundary, n_lag, t_lag, t_step, n_bootstraps; save_dir=false
    )    
    l_motif_an_timeseries = [[contribs[motif_class_num] for contribs in timeseries] for timeseries ∈ l_an_timeseries]
    peristimulus_results = map(test_sizes) do test_size
        peristimulus_start, peristimulus_stop = calculate_jitter_peristimulus_window(t_max_jitter, t_window, t_step, t_size ÷ 2)
        effs_and_sigs = map(1:n_resamples) do _
            l_timeseries_sample = rand(l_motif_an_timeseries, test_size)
            test_peristimulus_difference(l_timeseries_sample, peristimulus_start, peristimulus_stop)
        end
        mean_effect = mean(p.effect_size for p ∈ effs_and_sigs)
        proportion_rejected = mean(p.significance < α for p ∈ effs_and_sigs)
        std_effect = std(p.effect_size for p ∈ effs_and_sigs)
        (mean_effect=mean_effect, std_effect=std_effect, proportion_rejected=proportion_rejected, sample_size=test_size)
    end
    (l_an_timeseries, trialavg_raster, peristimulus_results)
end))

big_fig = Figure()

signal_raster = embedded_rand_motif(motif_class, n_size, t_size, -n_max_jitter:n_max_jitter, -t_max_jitter:t_max_jitter, n_max_jitter, t_max_jitter)
noise_raster = rand(size(signal_raster)...) .< noise_rate
raster = min.(signal_raster .+ noise_raster, 1)

f_signal = heatmap(signal_raster', axis=(xlabel="time", ylabel="neuron"))
f_noise = heatmap(noise_raster', axis=(xlabel="time", ylabel="neuron"))
f_raster = heatmap(raster', axis=(xlabel="time", ylabel="neuron"))
f_trialavg_raster = heatmap(trialavg_raster', axis=(xlabel="time", ylabel="neuron"))

results_df = DataFrame(peristimulus_results)
f_power = data(results_df) * mapping(:sample_size, :proportion_rejected) * visual(Scatter) |> draw

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

save(plotsdir(subdir,"power_vs_sample_size_$(motif_class).$(plot_ext)"), f_power)

(l_an_timeseries, trialavg_raster, peristimulus_results)
end
end

@save plotsdir(subdir, "results.jld2") prior_results_dict

