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
trials=100
n_bootstraps=20
n_resamples=50
n_test_points=20
α = 0.05
results_key = (; boundary=boundary, trials=trials, n_bootstraps=n_bootstraps, n_resamples=n_resamples, α=α, n_test_points=n_test_points)
subdir = if boundary isa Periodic
    "AN_$(trials)trials_$(n_bootstraps)bs_periodic_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
elseif boundary isa ZeroPadded
    "AN_$(trials)trials_$(n_bootstraps)bs_zeropad_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
else
    error("Unrecognized boundary condition for TriCorr")
end
mkpath(plotsdir(subdir))
save_all_trials_dir = plotsdir(subdir, "trials")
if save_all_trials_dir != false
    mkpath(save_all_trials_dir)
end

@warn "n_size set to 32"
let n_size = 32, t_size = 60,
    n_max_jitter = 3, t_max_jitter = 2,
    n_lag = 6, t_lag = 5, t_step=2,
    t_window = 2t_lag + 1,
    noise_rate = 0.2;

# # Middle p0
# n0_range = (n_max_jitter+1):(n_size-n_max_jitter); t0 = -t_max_jitter:t_max_jitter .+ (t_size ÷ 2)
# Latter half p0
n0_range = (n_max_jitter+1):(n_size-n_max_jitter)
t0_range = ((t_size ÷ 2) + t_max_jitter):(t_size - t_max_jitter)
@show n0_range t0_range

motif_class_range = 1:5

@threads for motif_class_num = motif_class_range
motif_class = offset_motif_numeral(motif_class_num)

test_sizes = 1:max(trials÷n_test_points,1):trials
get!(prior_results_dict, merge((motif_class=motif_class,), results_key), (begin
    l_an_timeseries, trialavg_raster = an_timeseries_across_jittered_trials(
        motif_class_num, n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter, trials, noise_rate, boundary, n_lag, t_lag, t_step, n_bootstraps; save_dir=save_all_trials_dir
    )    
    l_motif_an_timeseries = [[contribs[motif_class_num] for contribs in timeseries] for timeseries ∈ l_an_timeseries]
    peristimulus_results = map(test_sizes) do test_size
        peristimulus_start, peristimulus_stop = calculate_jitter_peristimulus_window(t0_range, t_max_jitter, t_window, t_step, t_size)
        effs_and_sigs = map(1:n_resamples) do _
            l_timeseries_sample = StatsBase.sample(l_motif_an_timeseries, test_size, replace=false)
            test_peristimulus_difference(l_timeseries_sample, peristimulus_start, peristimulus_stop)
        end
        mean_effect = mean(p.effect_size for p ∈ effs_and_sigs)
        proportion_rejected = mean(p.significance < α for p ∈ effs_and_sigs)
        std_effect = std(p.effect_size for p ∈ effs_and_sigs)
        (mean_effect=mean_effect, std_effect=std_effect, proportion_rejected=proportion_rejected, sample_size=test_size)
    end
    (l_an_timeseries, trialavg_raster, peristimulus_results)
end))
end

for motif_class_num = motif_class_range
    motif_class = offset_motif_numeral(motif_class_num)

    l_an_timeseries, trialavg_raster, peristimulus_results = 
    prior_results_dict[merge((motif_class=motif_class,), results_key)]

    signal_raster = embedded_rand_motif(motif_class, n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter)
    noise_raster = rand(size(signal_raster)...) .< noise_rate
    raster = min.(signal_raster .+ noise_raster, 1)

    f_signal = heatmap(signal_raster', axis=(xlabel="time", ylabel="neuron"))
    f_noise = heatmap(noise_raster', axis=(xlabel="time", ylabel="neuron"))
    f_raster = heatmap(raster', axis=(xlabel="time", ylabel="neuron"))
    f_trialavg_raster = heatmap(trialavg_raster', axis=(xlabel="time", ylabel="neuron"))

    results_df = DataFrame(peristimulus_results)
    plt_power = data(results_df) * mapping(:sample_size, :proportion_rejected) * visual(Scatter)
    f_power = draw(plt_power, axis=(; title=("Motif-class $(motif_class)")))

    f_motif_course = plot(1:t_step:(size(raster,2)-2t_lag), [a[motif_class_num] for a ∈ mean(l_an_timeseries)], axis=(xlabel="time", ylabel="avg motif $(motif_class) contrib"))
    f_motif_control = plot(1:t_step:(size(raster,2)-2t_lag), [a[14] for a ∈ mean(l_an_timeseries)], axis=(xlabel="time", ylabel="avg motif XIV contrib"))

    save(plotsdir(subdir,"signal_motif_$(motif_class)_AN.$(plot_ext)"), f_signal)
    save(plotsdir(subdir,"noise_motif_$(motif_class)_AN.$(plot_ext)"), f_noise)
    save(plotsdir(subdir,"raster_motif_$(motif_class)_AN.$(plot_ext)"), f_raster)
    save(plotsdir(subdir,"trialavg_raster_motif_$(motif_class)_AN.$(plot_ext)"), f_trialavg_raster)

    save(plotsdir(subdir,"motif_$(motif_class)_AN_timeseries_$(typeof(boundary)).$(plot_ext)"), f_motif_course)
    save(plotsdir(subdir,"motif_XIV_given_$(motif_class)_AN_timeseries_$(typeof(boundary)).$(plot_ext)"), f_motif_control)

    save(plotsdir(subdir,"power_vs_sample_size_$(motif_class).$(plot_ext)"), f_power)

end

@save plotsdir(subdir, "results.jld2") prior_results_dict n_size t_size n_max_jitter t_max_jitter n_lag t_lag t_step t_window noise_rate

end
