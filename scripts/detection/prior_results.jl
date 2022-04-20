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
using Colors

include(srcdir("roman_encode.jl"))
include(srcdir("peristimulus_testing.jl"))

function parse_results_key(; motif_class, boundary, trials, n_bootstraps, n_resamples, α, n_test_points)
    return (motif_class, boundary, trials, n_bootstraps, n_resamples, α, n_test_points)
end

@warn "n_size set to 32"
let n_size = 32, t_size = 60,
    n_max_jitter = 3, t_max_jitter = 2,
    n_lag = 6, t_lag = 5, t_step=2,
    t_window = 2t_lag + 1,
    noise_rate = 0.2,
    N_MOTIFS=14;

subdir = "multimotif_powers_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
mkpath(plotsdir(subdir))

all_saves = mapreduce(vcat, walkdir(projectdir("plots"))) do (root, dirs, files)
    basenames = filter(file -> occursin("jld2", file), files)
    joinpath.(Ref(root), basenames)
end
recent_save = sort!(all_saves; by=x -> stat(x).mtime)[end]

#@load projectdir("plots","AN_100trials_20bs_periodic_2022_04_15-114608", "results.jld2") prior_results_dict
@load recent_save prior_results_dict

keyed_motif_powers = map(keys(prior_results_dict) |> collect) do key
    signal_motif_class, boundary, trials, n_bootstraps, n_resamples, α, n_test_points = parse_results_key(; key...)

    test_sizes = 1:max(trials÷n_test_points,1):trials
    peristimulus_start, peristimulus_stop = calculate_jitter_peristimulus_window(t_max_jitter, t_window, t_step, t_size ÷ 2)

    l_trial_timeseries, trialavg_raster, peristimulus_results = prior_results_dict[key]
    l_trial_timeseries = reduce.(hcat, l_trial_timeseries)

    (key, mapreduce(vcat, axes(l_trial_timeseries[begin],1)) do target_motif_num
        mapreduce(vcat, test_sizes) do test_size
            effs_and_sigs = map(1:n_resamples) do _
                sampled_trials = StatsBase.sample(l_trial_timeseries, test_size, replace=false)
                target_motif_sampled_trials = map(sampled_trials) do trial
                    trial[target_motif_num, :]
                end
                test_peristimulus_difference(target_motif_sampled_trials, peristimulus_start, peristimulus_stop)
            end
            mean_effect = mean(p.effect_size for p ∈ effs_and_sigs)
            proportion_rejected = mean(p.significance < α for p ∈ effs_and_sigs)
            std_effect = std(p.effect_size for p ∈ effs_and_sigs)
            (mean_effect=mean_effect, std_effect=std_effect, proportion_rejected=proportion_rejected, sample_size=test_size, signal_motif=signal_motif_class, target_motif=offset_motif_numeral(target_motif_num))
        end
    end)
end

for (key, signaled_motif_powers) ∈ keyed_motif_powers
    signal_motif_class, boundary, trials, n_bootstraps, n_resamples, α, n_test_points = parse_results_key(; key...)

    powers_df = DataFrame(signaled_motif_powers)
    plt_powers = data(powers_df) * mapping(:sample_size, :proportion_rejected, color=:target_motif) * visual(Scatter)
    color_list = distinguishable_colors(N_MOTIFS, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    f_powers = draw(plt_powers, axis=(; title="Signal motif-class $(signal_motif_class)"), palettes=(color=color_list,))

    plt_effects = data(powers_df) * mapping(:sample_size, :mean_effect, color=:target_motif) * visual(Scatter)
    f_effects = draw(plt_effects, axis=(; title="Signal motif-class $(signal_motif_class)"), palettes=(color=color_list,))

    save(plotsdir(subdir,"power_detect_motif_$(signal_motif_class)_bd_$(boundary)_trials_$(trials)_bs_$(n_bootstraps)_rs_$(n_resamples)_α_$(α)_tp_$(n_test_points).$(plot_ext)"), f_powers)

    save(plotsdir(subdir,"effect_detect_motif_$(signal_motif_class)_bd_$(boundary)_trials_$(trials)_bs_$(n_bootstraps)_rs_$(n_resamples)_α_$(α)_tp_$(n_test_points).$(plot_ext)"), f_effects)


end

end #let


