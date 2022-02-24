
using Memoize

function test_peristimulus_difference(l_timeseries, start, stop)
    peristimulus = vcat([timeseries[start:stop] for timeseries ∈ l_timeseries]...)
    nonstimulus = vcat([[timeseries[begin:start-1]; timeseries[stop+1:end]] for timeseries ∈ l_timeseries]...)
    effect_size = mean(peristimulus) - mean(nonstimulus)
    significance = pvalue(MannWhitneyUTest(peristimulus, nonstimulus))
    return (effect_size=effect_size, significance=significance)
end

function get_motif_t_bounds(motif_class)
    triplet = TripleCorrelations.motif_examples[motif_class]
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

function make_an_timeseries(raster, boundary, n_lag, t_lag, t_step; t_window=2t_lag+1, n_bootstraps)
    N,T = size(raster)
    map(1:t_step:(T-t_window)) do t_start
        bootstrap_normed_sequence_classes(raster[:,t_start:t_start+t_window], boundary, n_lag, t_lag; 
            n_bootstraps=n_bootstraps
        )
    end
end

function detect_an_across_trials(signal_raster, trials, noise_rate, boundary, n_lag, t_lag, t_step, n_bootstraps)

    @showprogress map(1:trials) do _
        noise_raster = rand(size(signal_raster)...) .< noise_rate
        raster = (signal_raster .+ noise_raster) .> 0
        make_an_timeseries(raster, boundary, n_lag, t_lag, t_step; n_bootstraps=n_bootstraps)
    end

end