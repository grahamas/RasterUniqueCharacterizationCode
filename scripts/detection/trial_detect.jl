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
include(srcdir("filename_utils.jl"))
include(srcdir("plotting.jl"))

trials_epoch_tricorrs = let conditioned_on = Constituents(),
    postproc! = demean!;

#n_signals = 1
norming="rate_divide"
boundary = PeriodicExtended(50)
boundary_width = get_boundary_width(boundary)
#n_trials=35
#n_tests=25
α = 0.05 / 14

MEA_n = 100
MEA_t = 101
#MEA_spikes = 261

raster_size = (MEA_n, (MEA_t*2)+(boundary_width*4))
#signal_jitters = (2, 2)
#lag_extents = (5, 5)
noise_rate = MEA_spikes / (MEA_n * MEA_t)

only_rate_raster = zeros(Bool, MEA_n, MEA_t + boundary_width*2)
raster_meat = slice_meat(only_rate_raster, boundary)
raster_meat[1:MEA_spikes] .= true
if boundary_width != 0
    expected_spikes = round(Int, MEA_spikes * (2 * boundary_width) / (MEA_t))
    half_spikes = floor(Int, MEA_spikes / 2)
    flank1, flank2 = slice_flanks(only_rate_raster, boundary)
    flank1[1:half_spikes] .= true
    flank2[1:(expected_spikes - half_spikes)] .= true
end

@show count(only_rate_raster)
@show count(slice_meat(only_rate_raster, boundary))

only_rate_raster[1:MEA_spikes] .= true
precalced_postproc! = precalculate(postproc!, IndBernoulli(), conditioned_on, [only_rate_raster], boundary, lag_extents)

outputdir = plotsdir("particular_trial_detects_constituents_$(floor(Int, time() * 1000))")
mkpath(outputdir)


### Trial detection

trial_results = @time run_peristimulus_tests(
    raster_size, boundary, lag_extents, 
    n_signals, signal_jitters, noise_rate,
    n_trials, n_tests,
    precalced_postproc!, α
)

trial_results.signal_motif_numeral = offset_motif_numeral.(trial_results.signal_motif)
trial_results.detected_motif_numeral = offset_motif_numeral.(trial_results.detected_motif)
plt = data(trial_results) * mapping(:signal_motif_numeral => sorter(1:14 .|> offset_motif_numeral) => "signal_motif", :detected_motif_numeral => sorter(1:14 .|> offset_motif_numeral) => "detected motif", color=:proportion_detected => "prop. tests positive") * visual(markersize=29)

axis = (; title = "$(n_trials) trials | $(n_tests) tests | $(n_signals) signals")
fig = draw(plt; axis)

save(joinpath(outputdir, "crossmotif_significances_$(fn2str(postproc!))_$(obj2str(conditioned_on))_trials$(n_trials)_tests$(n_tests)_signals$(n_signals)_spikes$(MEA_spikes)_size$(raster_size)_lags$(lag_extents).png"), fig)


# Example raster

motif_class = "VII"
raster_size = (MEA_n, (MEA_t*2)+(boundary_width*4))
n_size, t_size = raster_size
n_max_jitter, t_max_jitter = signal_jitters
boundary_width = get_boundary_width(boundary)

n0_range = (n_max_jitter+1):(n_size-n_max_jitter)
t0_range = ((t_size ÷ 2) + t_max_jitter+boundary_width):(t_size - t_max_jitter-boundary_width)
base_node_ranges = (n0_range, t0_range)

epochs = [1:t0_range[begin]-t_max_jitter-boundary_width,(t0_range[begin]-t_max_jitter-boundary_width+1):raster_size[end]]

raster = embedded_rand_motif(motif_class, raster_size, base_node_ranges, signal_jitters)
for _ in 1:(n_signals-1)
    raster .|= embedded_rand_motif(motif_class, raster_size, base_node_ranges, signal_jitters)
end
signal_raster = copy(raster)
for epoch in epochs
    epoch_raster = view_slice_last(raster, epoch)
    noise_spikes = round(Int, prod(size(epoch_raster)) * noise_rate)
    spikes_left = noise_spikes - count(epoch_raster) 
    while spikes_left > 0
        epoch_noise_raster = fixed_noise_raster(epoch_raster, noise_rate, boundary)
        epoch_raster .|= epoch_noise_raster
        spikes_left = noise_spikes - count(epoch_raster)
    end
end

raster_resolution = (800, 600)
signal_fig = Figure(resolution=raster_resolution)
signal_axis = Axis(signal_fig[1,1])
plot_raster!(signal_axis, signal_raster)

raster_fig = Figure(resolution=raster_resolution)
raster_axis = Axis(raster_fig[1,1])
plot_raster!(raster_axis, raster)

save(joinpath(outputdir, "signals_signals$(n_signals)_spikes$(MEA_spikes)_size$(raster_size)_lags$(lag_extents).png"), signal_fig)
save(joinpath(outputdir, "embedded_signals$(n_signals)_spikes$(MEA_spikes)_size$(raster_size)_lags$(lag_extents).png"), raster_fig)


####### Bar and whisker plot

@show n_trials
trials_epoch_tricorrs = map(1:n_trials) do _
    raster = embedded_rand_motif(motif_class, raster_size, base_node_ranges, signal_jitters)
    for _ in 1:(n_signals-1)
        raster .|= embedded_rand_motif(motif_class, raster_size, base_node_ranges, signal_jitters)
    end  
    for epoch in epochs
        epoch_raster = view_slice_last(raster, epoch)
        set_noise_raster!(epoch_raster, boundary, noise_rate)
    end
    epoch_tricorrs = calculate_trial_epochs(raster, boundary, lag_extents, epochs, precalced_postproc!)
    epoch_tricorrs
end

nts = mapreduce(vcat, enumerate(trials_epoch_tricorrs)) do (i_trial, trial)
    control_epochs = trial[:,1]
    signal_epochs = trial[:,2]
    motifs = offset_motif_numeral.(1:14)
    mapreduce(vcat, zip(control_epochs, signal_epochs, motifs)) do (c, s, m)
        [
            (value=c, epoch="control", trial=i_trial, motif=m),
            (value=s, epoch="signal", trial=i_trial, motif=m)
        ]
    end
end
df = DataFrame(nts)

epoch_plot = data(df) * mapping(:motif, :value, dodge=:epoch, color=:epoch) * visual(BoxPlot)
epoch_fig = draw(epoch_plot, axis=(title="Detecting motif $(motif_class)", ylabel="A-T"))

save(joinpath(outputdir, "epoch_detection_constituents_motif_$(motif_class)_signals$(n_signals)_spikes$(MEA_spikes)_size$(raster_size)_lags$(lag_extents).png"), epoch_fig)


#############

trials_epoch_tricorrs
end
