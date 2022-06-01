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

if !@isdefined(prior_results_dict) || force_redef
    prior_results_dict = Dict()
end
force_redef = false

contribution_fn_dict = Dict(
	"constituent_divide" => sequence_classes_divide_E_given_constituents,
	"rate_divide" => sequence_classes_divide_E_given_rate,
    "actual" => sequence_class_tricorr
)

n_signals = 10
norming="rate_divide"
N_MOTIFS=14
boundary = Periodic()#Extended(50)
boundary_width = boundary isa PeriodicExtended ? boundary.boundary : 0
trials=15
n_resamples=2
n_test_points=2
α = 0.05 / 14

MEA_n = 100
MEA_t = 101
MEA_spikes = 167
results_key = (; boundary=boundary, trials=trials, n_resamples=n_resamples, α=α, n_test_points=n_test_points)
subdir = if boundary isa Periodic
    "$(norming)_$(trials)trials_$(n_signals)sigs_MEA_$(MEA_n)x$(MEA_t)_$(MEA_spikes)spk_periodic_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
elseif boundary isa ZeroPadded
    "$(norming)_$(trials)trials_$(n_signals)sigs_MEA_$(MEA_n)x$(MEA_t)_$(MEA_spikes)spk_zeropad_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
elseif boundary isa PeriodicExtended
    "$(norming)_$(trials)trials_$(n_signals)sigs_MEA_$(MEA_n)x$(MEA_t)_$(MEA_spikes)spk_PeriodicExtended_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
else
    error("Unrecognized boundary condition for TriCorr")
end
mkpath(plotsdir(subdir))
save_all_trials_dir = false#plotsdir(subdir, "trials")
if save_all_trials_dir != false
    mkpath(save_all_trials_dir)
end

#n_size = 32; t_size = 60; n_max_jitter = 3; t_max_jitter = 2; n_lag = 12; t_lag = 10; t_step=2; noise_rate = 0.3; motif_class_num = 1; save_all_trials_dir=false

@warn "n_size set to 32"
let n_size = MEA_n, t_size = (MEA_t*2)+(boundary_width*4),
    n_max_jitter = 2, t_max_jitter = 3,
    n_lag = 10, t_lag = 8, t_step=2,
    noise_rate = MEA_spikes / (MEA_n * MEA_t);

contribution_fn = contribution_fn_dict[norming]

# # Middle p0
# n0_range = (n_max_jitter+1):(n_size-n_max_jitter); t0 = -t_max_jitter:t_max_jitter .+ (t_size ÷ 2)
# Latter half p0
n0_range = (n_max_jitter+1):(n_size-n_max_jitter)
t0_range = ((t_size ÷ 2) + t_max_jitter+boundary_width):(t_size - t_max_jitter-boundary_width)

motif_class_range = 1:14

@threads for motif_class_num = motif_class_range
motif_class = offset_motif_numeral(motif_class_num)

# trials_epoch_tricorrs:
# [trial : [epochs : motif x epoch]]]

test_sizes = max(trials÷n_test_points,15):trials÷n_test_points:trials
@assert length(test_sizes) > 0
@show merge((motif_class=motif_class,), results_key)
get!(prior_results_dict, merge((motif_class=motif_class,), results_key), (begin
    @info "Motif class $(motif_class) signal..."
    trials_epoch_tricorrs, trialavg_raster = @time jittered_trials_epochs(
        motif_class_num, n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter, trials, noise_rate, n_signals, boundary, (n_lag, t_lag), contribution_fn; save_dir=save_all_trials_dir
    )
    @info "done ($motif_class). Statistics..."
    peristimulus_results = @showprogress map(test_sizes) do test_size
        effs_and_sigs = mapreduce(vcat, 1:n_resamples) do _
            these_trials_epoch_tricorrs, _ = jittered_trials_epochs(
                motif_class_num, n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter, test_size, noise_rate, n_signals, boundary, (n_lag, t_lag), contribution_fn; save_dir=save_all_trials_dir
            )
            test_epoch_difference(these_trials_epoch_tricorrs)
        end
        mean_effect = mean(p.effect_size for p ∈ effs_and_sigs)
        proportion_rejected = mean(p.significance .< α for p ∈ effs_and_sigs)
        #std_effect = std(p.effect_size for p ∈ effs_and_sigs)
        (mean_effect=mean_effect, #std_effect=std_effect, 
        proportion_rejected=proportion_rejected, sample_size=test_size)
    end
    # peristimulus results are each motif contrib vectors
    (trials_epoch_tricorrs, trialavg_raster, peristimulus_results)
end))
end # throads for

@show "Made prior results..."

for motif_class_num = motif_class_range
    motif_class = offset_motif_numeral(motif_class_num)
    @show "Plots for Motif Class $motif_class"

    l_trials_epochs, trialavg_raster, peristimulus_results = 
    prior_results_dict[merge((motif_class=motif_class,), results_key)]

    # Example rasters (signal, noise, combined, trial average)
    signal_raster = embedded_rand_motif(motif_class, n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter)
	for _ in 1:(n_signals-1)
		signal_raster .|= embedded_rand_motif(motif_class, n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter)
	end
    noise_raster = rand(size(signal_raster)...) .< noise_rate
    raster = min.(signal_raster .+ noise_raster, 1)

    f_target = Figure()
    ax = Axis(f_target[1,1], xlabel="epoch", ylabel="A/E Motif $(motif_class)")
    lines!.(Ref(ax), Ref([1, 2]), [l_trials_epochs[i][motif_class_num,:] for i in 1:size(l_trials_epochs,1)])
    epoch1 = [l_trials_epochs[i][motif_class_num,1] for i ∈ 1:size(l_trials_epochs,1)]
    epoch2 = [l_trials_epochs[i][motif_class_num,2] for i ∈ 1:size(l_trials_epochs,1)]
    boxplot!([[1 for _ ∈ 1:length(epoch1)]..., [2 for _ ∈ 1:length(epoch2)]...], [epoch1..., epoch2...])
    ax.xticks = ([1,2], ["noise","signal"])

    f_signal = heatmap(signal_raster', axis=(xlabel="time", ylabel="neuron"))
    f_noise = heatmap(noise_raster', axis=(xlabel="time", ylabel="neuron"))
    f_raster = heatmap(raster', axis=(xlabel="time", ylabel="neuron"))
    f_trialavg_raster = heatmap(trialavg_raster', axis=(xlabel="time", ylabel="neuron"))

    epochs_df = DataFrame(mapreduce(vcat, l_trials_epochs) do trial
        @assert size(trial) == (14,2)
        mapreduce(vcat, axes(trial, 1)) do i
            [(
                signal_motif=motif_class, 
                detect_motif=offset_motif_numeral(i),
                condition=:nonstim,
                contrib=trial[i,1]
            ), 
            (
                signal_motif=motif_class, 
                detect_motif=offset_motif_numeral(i),
                condition=:stim,
                contrib=trial[i,2]
            )]
        end
    end)
    plt_epochs = data(epochs_df) * mapping(:contrib, color=:condition, layout=:detect_motif) * histogram(; bins=15)
    f_epochs = draw(plt_epochs)#, axis=(; title=("Motif-class $(motif_class) signal")))

    # Power timeseries
    peristimulus_results_by_motif = DataFrame(mapreduce(vcat, peristimulus_results) do result
        map(axes(result.mean_effect,1)) do i_motif
            (
                signal_motif=motif_class,
                detect_motif=TripleCorrelations.offset_motif_numeral(i_motif),
                mean_effect=result.mean_effect[i_motif], 
                #std_effect=result.std_effect[i_motif],
                proportion_rejected=result.proportion_rejected[i_motif],
                sample_size=result.sample_size
            )
        end
    end)
    plt_power = data(peristimulus_results_by_motif) * mapping(:sample_size, :proportion_rejected, color=:detect_motif) * (visual(Scatter, markersize=0.1, strokewidth=0) + smooth(span=0.9, degree=2))
    color_list = distinguishable_colors(N_MOTIFS, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    f_power = draw(plt_power, palettes=(color=color_list,))#, axis=(; title="Signal motif-class $(motif_class)"))

    plt_effect = data(peristimulus_results_by_motif) * mapping(:sample_size, :mean_effect, color=:detect_motif) * (visual(Scatter, markersize=0.1, strokewidth=0) + smooth(span=0.9, degree=2))
    color_list = distinguishable_colors(N_MOTIFS, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    f_effect = draw(plt_effect, palettes=(color=color_list,))#, axis=(; title="Signal motif-class $(motif_class)"))

    # f_motif_course = plot(1:2, [a[motif_class_num] for a ∈ mean(l_trials_epochs)], axis=(xlabel="time", ylabel="avg motif $(motif_class) contrib"))
    # f_motif_control = plot(1:2, [a[14] for a ∈ mean(l_trials_epochs)], axis=(xlabel="time", ylabel="avg motif XIV contrib"))

    # save(plotsdir(subdir,"motif_$(motif_class)_AN_timeseries_$(typeof(boundary)).$(plot_ext)"), f_motif_course)
    # save(plotsdir(subdir,"motif_XIV_given_$(motif_class)_AN_timeseries_$(typeof(boundary)).$(plot_ext)"), f_motif_control)

    save(plotsdir(subdir,"target_$(motif_class)_epoch_AE.$(plot_ext)"), f_target)

    save(plotsdir(subdir,"signal_motif_$(motif_class)_AN.$(plot_ext)"), f_signal)
    save(plotsdir(subdir,"noise_motif_$(motif_class)_AN.$(plot_ext)"), f_noise)
    save(plotsdir(subdir,"raster_motif_$(motif_class)_AN.$(plot_ext)"), f_raster)
    save(plotsdir(subdir,"trialavg_raster_motif_$(motif_class)_AN.$(plot_ext)"), f_trialavg_raster)

    save(plotsdir(subdir,"power_vs_sample_size_$(motif_class).$(plot_ext)"), f_power)
    save(plotsdir(subdir,"effect_vs_sample_size_$(motif_class).$(plot_ext)"), f_effect)

    save(plotsdir(subdir, "motif_contrib_distributions_$(motif_class).$(plot_ext)"), f_epochs)

end

@save plotsdir(subdir, "results.jld2") prior_results_dict n_size t_size n_max_jitter t_max_jitter n_lag t_lag t_step noise_rate n_signals

end
