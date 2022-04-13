
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

function calculate_jitter_peristimulus_window(t_jitter, t_window, t_step, motif_t_target)
    window_start_dxs = 1:t_step:motif_t_target
    window_stop_dxs = (1+t_window):t_step:motif_t_target+t_jitter
    first_touch = findfirst(window_stop_dxs .>= (motif_t_target - t_jitter)) # Must have entire motif
    last_touch = findlast(window_start_dxs .<= (motif_t_target + t_jitter))
    return (first_touch, last_touch)
end

function make_an_timeseries(raster, boundary, n_lag, t_lag, t_step; t_window=2t_lag+1, n_bootstraps)
    N,T = size(raster)
    map(1:t_step:(T-t_window)) do t_start
        bootstrap_normed_sequence_classes(raster[:,t_start:t_start+t_window], boundary, n_lag, t_lag; 
            n_bootstraps=n_bootstraps
        )
    end
end

function detect_an_across_trials(motif_class_num::Int, signal_raster::Array, trials, noise_rate, boundary, n_lag, t_lag, t_step, n_bootstraps; save_dir=false)
    motif_class = roman_encode_zero(motif_class_num)
    trialavg_raster = zeros(Float64, size(signal_raster)...)
    l_an_timeseries = @showprogress map(1:trials) do trial_num
        noise_raster = rand(size(signal_raster)...) .< noise_rate
        raster = Array{Bool}((signal_raster .+ noise_raster) .> 0)
        trialavg_raster += raster
        an_timeseries = make_an_timeseries(raster, boundary, n_lag, t_lag, t_step; n_bootstraps=n_bootstraps)
        if save_dir != false
            f_signal = heatmap(signal_raster', axis=(xlabel="time", ylabel="neuron"))
            f_noise = heatmap(noise_raster', axis=(xlabel="time", ylabel="neuron"))
            f_raster = heatmap(raster', axis=(xlabel="time", ylabel="neuron"))
            f_motif_course = plot([a[motif_class_num] for a ∈ an_timeseries], axis=(xlabel="time", ylabel="avg motif $(motif_class) contrib"))
            save(joinpath(save_dir ,"signal_motif_$(motif_class)_AN_$(trial_num).$(plot_ext)"), f_signal)
            save(joinpath(save_dir ,"noise_motif_$(motif_class)_AN_$(trial_num).$(plot_ext)"), f_noise)
            save(joinpath(save_dir ,"raster_motif_$(motif_class)_AN_$(trial_num).$(plot_ext)"), f_raster)
            save(joinpath(save_dir ,"motif_$(motif_class)_AN_$(trial_num)_timeseries_$(typeof(boundary)).$(plot_ext)"), f_motif_course)
        end
        an_timeseries
    end
    trialavg_raster ./= trials
    return (l_an_timeseries, trialavg_raster)
end

function detect_feedback_tricorr(tricorr::TripleCorrelation)

end


function embedded_rand_motif(motif_class, n_size, t_size, n0_range::AbstractArray, t0_range::AbstractArray, n_max_jitter, t_max_jitter)
    raster = zeros(Bool, n_size, t_size)
    motif_coords = TripleCorrelations.rand_motif(motif_class, n0_range, t0_range, n_max_jitter, t_max_jitter)
    for motif_coord in motif_coords
        raster[CartesianIndex((motif_coord .+ (n_size ÷ 2, t_size ÷ 2)) .% (n_size, t_size) .+ (1, 1))] = 1
    end
    return raster
end

@memoize function detect_an_across_jittered_trials(motif_class_num::Int, n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter, trials, noise_rate, boundary, n_lag, t_lag, t_step, n_bootstraps; save_dir=false)
    motif_class = offset_motif_numeral(motif_class_num)
    trialavg_raster = zeros(Float64, n_size, t_size)
    l_an_timeseries = @showprogress map(1:trials) do trial_num
        signal_raster = embedded_rand_motif(motif_class, n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter)
        noise_raster = rand(size(signal_raster)...) .< noise_rate
        raster = Array{Bool}((signal_raster .+ noise_raster) .> 0)
        trialavg_raster += raster
        an_timeseries = make_an_timeseries(raster, boundary, n_lag, t_lag, t_step; n_bootstraps=n_bootstraps)
        if save_dir != false
            f_signal = heatmap(signal_raster', axis=(xlabel="time", ylabel="neuron"))
            f_noise = heatmap(noise_raster', axis=(xlabel="time", ylabel="neuron"))
            f_raster = heatmap(raster', axis=(xlabel="time", ylabel="neuron"))
            f_motif_course = plot(1:t_step:(size(raster,2)-2t_lag), [a[motif_class_num] for a ∈ an_timeseries], axis=(xlabel="time", ylabel="avg motif $(motif_class) contrib"))
            save(joinpath(save_dir ,"signal_motif_$(motif_class)_AN_$(trial_num).$(plot_ext)"), f_signal)
            save(joinpath(save_dir ,"noise_motif_$(motif_class)_AN_$(trial_num).$(plot_ext)"), f_noise)
            save(joinpath(save_dir ,"raster_motif_$(motif_class)_AN_$(trial_num).$(plot_ext)"), f_raster)
            save(joinpath(save_dir ,"motif_$(motif_class)_AN_$(trial_num)_timeseries_$(typeof(boundary)).$(plot_ext)"), f_motif_course)
        end
        an_timeseries
    end
    trialavg_raster ./= trials
    return (l_an_timeseries, trialavg_raster)
end