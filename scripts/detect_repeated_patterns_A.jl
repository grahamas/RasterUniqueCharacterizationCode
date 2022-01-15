@quickactivate "RasterUniqueCharacterizationCode"

using TripleCorrelations
using HypothesisTests, Random, Statistics
using ProgressMeter
using GLMakie

#defines motif_examples dict
include("../../TripleCorrelations/test/data/motif_examples.jl")
include("../../TripleCorrelations/test/src/helpers.jl")
include("../src/roman_encode.jl")

if !@isdefined(a_timeseries_dict)
    a_timeseries_dict = Dict()
end

if !@isdefined(peristimulus_a_results_dict)
    peristimulus_a_results_dict = Dict()
end


function make_a_timeseries(raster, n_lag, t_lag, t_step; t_window=2t_lag+1)
    N,T = size(raster)
    map(1:t_step:(T-t_window)) do t_start
        sequence_class_tricorr(raster[:,t_start:t_start+t_window], n_lag, t_lag)
    end
end

function test_peristimulus_difference(l_timeseries, motif_class_num, start, stop)
    peristimulus = vcat([[contribs[motif_class_num] for contribs ∈ timeseries[start:stop]] for timeseries ∈ l_timeseries]...)
    nonstimulus = vcat([[contribs[motif_class_num] for contribs ∈ [timeseries[begin:start-1]; timeseries[stop+1:end]]] for timeseries ∈ l_timeseries]...)
    effect_size = mean(peristimulus) - mean(nonstimulus)
    significance = pvalue(MannWhitneyUTest(peristimulus, nonstimulus))
    return (N=length(l_timeseries), effect_size=effect_size, significance=significance)
end

function get_motif_t_bounds(motif_class)
    triplet = motif_examples[motif_class]
    time_collapsed = dropdims(sum(triplet, dims=1), dims=1)
    (findfirst(time_collapsed .> 0), findlast(time_collapsed .> 0))
end


function calculate_peristimulus_window(t_pad, t_window, t_step, motif_t_start, motif_t_stop, motif_total_width=3)
    total_length = 2t_pad + motif_total_width
    start_dxs = 1:t_step:(total_length - t_pad)
    stop_dxs = (1+t_window):t_step:total_length
    first_touch = findfirst(stop_dxs .>= (t_pad + motif_t_stop)) # Must have entire motif
    last_touch = findlast(start_dxs .<= (t_pad + motif_t_start))
    return (first_touch+1, last_touch-1)
end






motif_class_num = 5
motif_class = roman_encode(motif_class_num)
a_timeseries_dict[motif_class], peristimulus_a_results_dict[motif_class] = let n_pad = 10, t_pad = 30, 
    n_reps = 1, t_reps = 1,
    n_lag = 6, t_lag = 5,
    snippets=500, t_step=2,
    t_window = 2t_lag + 1,
    noise_rate = 0.2, noise_calcs=2;

signal_raster = repeat_padded_motif(motif_class, n_pad, t_pad, n_reps, t_reps)
l_a_timeseries = @showprogress map(1:snippets) do i_snippet
    noise_raster = rand(size(signal_raster)...) .< noise_rate
    raster = min.(signal_raster .+ noise_raster, 1)
    make_a_timeseries(raster, n_lag, t_lag, t_step)
end
test_sizes = 1:max(snippets÷10,1):snippets
peristimulus_results = map(test_sizes) do test_size
    l_timeseries_sample = rand(l_a_timeseries, test_size)
    motif_t_bounds = get_motif_t_bounds(motif_class)
    peristimulus_start, peristimulus_stop = calculate_peristimulus_window(t_pad, t_window, t_step, motif_t_bounds...)
    test_peristimulus_difference(l_timeseries_sample, motif_class_num, peristimulus_start, peristimulus_stop)
end

signal_raster = repeat_padded_motif(motif_class, n_pad, t_pad, n_reps, t_reps)
noise_raster = rand(size(signal_raster)...) .< noise_rate
raster = min.(signal_raster .+ noise_raster, 1)

f_signal = heatmap(signal_raster')
f_noise = heatmap(noise_raster')
f_raster = heatmap(raster')

f_motif_course = plot([a[motif_class_num] for a ∈ mean(l_a_timeseries)])
f_motif_control = plot([a[14] for a ∈ mean(l_a_timeseries)])

save(plotsdir("signal_motif_A_$(motif_class).png"), f_signal)
save(plotsdir("noise_motif_A_$(motif_class).png"), f_noise)
save(plotsdir("raster_motif_A_$(motif_class).png"), f_raster)
save(plotsdir("motif_A_timeseries_$(motif_class).png"), f_motif_course)
save(plotsdir("motif_XIV_A_timeseries_given_$(motif_class).png"), f_motif_control)

(l_a_timeseries, peristimulus_results)
end

