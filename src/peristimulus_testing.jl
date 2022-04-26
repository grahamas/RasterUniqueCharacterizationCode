
using Memoize

function test_peristimulus_difference(l_timeseries::AbstractVector{<:AbstractVector}, start, stop)
    peristimulus = vcat([timeseries[start:stop] for timeseries ∈ l_timeseries]...)
    nonstimulus = vcat([[timeseries[begin:start-1]; timeseries[stop+1:end]] for timeseries ∈ l_timeseries]...)
    effect_size = mean(peristimulus) - mean(nonstimulus)
    significance = pvalue(UnequalVarianceTTest(peristimulus, nonstimulus))
    return (effect_size=effect_size, significance=significance)
end

function test_epoch_difference(l_trials::AbstractVector{<:AbstractVector})
    stimulus = [trial[2] for trial ∈ l_trials]
    nonstimulus = [trial[1] for trial ∈ l_trials]
    effect_size = mean(stimulus) - mean(nonstimulus)
    significance = pvalue(UnequalVarianceTTest(stimulus, nonstimulus))
    return (effect_size=effect_size, significance=significance)
end

function test_epoch_difference(l_trials::AbstractVector{<:AbstractMatrix})
    # matrix of contributions
    n_motifs = size(first(l_trials),1)
    @assert n_motifs == N_MOTIFS
    hypothesis_pairs = map(1:n_motifs) do i_motif
        stimulus = [trial[i_motif,2] for trial ∈ l_trials]
        nonstimulus = [trial[i_motif,1] for trial in l_trials]
        (stimulus, nonstimulus)
    end
    effect_size = hypothesis_pairs .|> (x) -> mean(x[1]) - mean(x[2])
    significance = [pvalue(UnequalVarianceTTest(pair...)) for pair in hypothesis_pairs]
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

function calculate_jitter_peristimulus_window(t0_range, t_jitter, t_window, t_step, t_size)
    motif_start = t0_range[begin]
    motif_stop = t0_range[end]
    window_start_dxs = 1:t_step:(t_size-t_window)
    window_stop_dxs = (1+t_window):t_step:t_size
    first_touch = findfirst(window_stop_dxs .>= (motif_start - t_jitter)) # Must have entire motif
    last_touch = findlast(window_start_dxs .<= (motif_stop + t_jitter))
    return (first_touch, last_touch)
end

# unction calculate_jitter_peristimulus_window(t_jitter, t_window, t_step, motif_t_target)
#     window_start_dxs = 1:t_step:motif_t_target
#     window_stop_dxs = (1+t_window):t_step:motif_t_target+t_jitter
#     first_touch = findfirst(window_stop_dxs .>= (motif_t_target - t_jitter)) # Must have entire motif
#     last_touch = findlast(window_start_dxs .<= (motif_t_target + t_jitter))
#     return (first_touch, last_touch)
# end

function make_an_timeseries(raster, boundary, lag_extents, t_step; t_window=2t_lag+1, n_bootstraps)
    N,T = size(raster)
    map(1:t_step:(T-t_window)) do t_start
        bootstrap_normed_sequence_classes(raster[:,t_start:t_start+t_window], boundary, (lag_extents); 
            n_bootstraps=n_bootstraps
        )
    end
end

function detect_an_across_trials(motif_class_num::Int, signal_raster::Array, trials, noise_rate, boundary, lag_extents, t_step, n_bootstraps; save_dir=false)
    motif_class = offset_motif_numeral(motif_class_num)
    trialavg_raster = zeros(Float64, size(signal_raster)...)
    l_an_timeseries = @showprogress map(1:trials) do trial_num
        noise_raster = rand(size(signal_raster)...) .< noise_rate
        raster = Array{Bool}((signal_raster .+ noise_raster) .> 0)
        trialavg_raster += raster
        an_timeseries = make_an_timeseries(raster, boundary, lag_extents, t_step; n_bootstraps=n_bootstraps)
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


function embedded_rand_motif(motif_class, n_size, t_size, n0_range::AbstractArray, t0_range::AbstractArray, n_max_jitter, t_max_jitter)
    raster = zeros(Bool, n_size, t_size)
    motif_coords = TripleCorrelations.rand_motif(motif_class, n0_range, t0_range, n_max_jitter, t_max_jitter)
    for motif_coord in motif_coords
        raster[CartesianIndex(motif_coord)] = 1
    end
    return raster
end

function an_timeseries_across_jittered_trials(motif_class_num::Int, n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter, trials, noise_rate, boundary, lag_extents, t_step, n_bootstraps; save_dir=false)
    motif_class = offset_motif_numeral(motif_class_num)
    trialavg_raster = zeros(Float64, n_size, t_size)
    l_an_timeseries = @showprogress map(1:trials) do trial_num
        signal_raster = embedded_rand_motif(motif_class, n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter)
        noise_raster = rand(size(signal_raster)...) .< noise_rate
        raster = Array{Bool}((signal_raster .+ noise_raster) .> 0)
        trialavg_raster += raster
        an_timeseries = make_an_timeseries(raster, boundary, lag_extents, t_step; n_bootstraps=n_bootstraps)
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

function calculate_trial_epochs(raster, boundary, lag_extents, epochs; n_bootstraps)
    mapreduce(hcat, epochs) do epoch
        @show epoch size(raster[:,epoch]) sum(raster[:,epoch])
        bootstrap_normed_sequence_classes(raster[:,epoch], boundary, (lag_extents); n_bootstraps=n_bootstraps)
    end
end

@inline function view_slice_last(arr::AbstractArray{T,N}, dx) where {T,N}
    view(arr, ntuple(_ -> Colon(), N - 1)..., dx)
end

function fixed_noise_raster(dims, noise_rate, boundary::PeriodicExtended)
    n_bdry_ones = round(Int, noise_rate * prod(dims[1:end-1]) * boundary.boundary)
    n_meat_ones = round(Int, noise_rate * prod(dims[1:end-1]) * (dims[end] - 2(boundary.boundary)))
    noise_raster = zeros(Bool, dims...)
    bdry_begin = view_slice_last(noise_raster, 1:boundary.boundary)
    meat = view_slice_last(noise_raster, boundary.boundary+1:(dims[end]-boundary.boundary))
    bdry_end = view_slice_last(noise_raster, (dims[end]-boundary.boundary+1):dims[end])
    bdry_begin[1:n_bdry_ones] .= 1
    meat[1:n_meat_ones] .= 1
    bdry_end[1:n_bdry_ones] .= 1
    shuffle!(bdry_begin)
    shuffle!(meat)
    shuffle!(bdry_end)
    return noise_raster
end

function fixed_noise_raster(dims, noise_rate, boundary)
    n_ones = round(Int, prod(dims) * noise_rate)
    noise_raster = zeros(Bool, dims...)
    noise_raster[1:n_ones] .= 1
    shuffle!(noise_raster)
    return noise_raster
end

function jittered_trials_epochs(motif_class_num::Int, 
        n_size, t_size, 
        n0_range, t0_range, 
        n_max_jitter, t_max_jitter, 
        trials, noise_rate, boundary, 
        lag_extents, n_bootstraps; save_dir=false
    )
    motif_class = offset_motif_numeral(motif_class_num)
    trialavg_raster = zeros(Float64, n_size, t_size)
    trials_epoch_tricorrs = map(1:trials) do trial_num
        raster = embedded_rand_motif(motif_class, n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter)
        noise_raster = fixed_noise_raster(size(raster), noise_rate, boundary)
        raster .|= noise_raster
        @show raster
        trialavg_raster += raster
        @assert t0_range[begin] > 1+t_max_jitter
        epochs = [1:t0_range[begin]-t_max_jitter,(t0_range[begin]-t_max_jitter+1):size(raster)[end]]
        @show length.(epochs)
        epoch_tricorrs = calculate_trial_epochs(raster, boundary, lag_extents, epochs; n_bootstraps=n_bootstraps)
        if save_dir != false
            f_signal = heatmap(signal_raster', axis=(xlabel="time", ylabel="neuron"))
            f_noise = heatmap(noise_raster', axis=(xlabel="time", ylabel="neuron"))
            f_raster = heatmap(raster', axis=(xlabel="time", ylabel="neuron"))
            f_motif_course = plot(1:length(epochs), [a[motif_class_num] for a ∈ epoch_tricorrs], axis=(xlabel="time", ylabel="avg motif $(motif_class) contrib"))
            save(joinpath(save_dir ,"signal_motif_$(motif_class)_AN_$(trial_num).$(plot_ext)"), f_signal)
            save(joinpath(save_dir ,"noise_motif_$(motif_class)_AN_$(trial_num).$(plot_ext)"), f_noise)
            save(joinpath(save_dir ,"raster_motif_$(motif_class)_AN_$(trial_num).$(plot_ext)"), f_raster)
            save(joinpath(save_dir ,"motif_$(motif_class)_AN_$(trial_num)_timeseries_$(typeof(boundary)).$(plot_ext)"), f_motif_course)
        end
        epoch_tricorrs
    end
    trialavg_raster ./= trials
    return (trials_epoch_tricorrs, trialavg_raster)
end