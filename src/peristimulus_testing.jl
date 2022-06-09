
using Memoize, ThreadsX

get_boundary_width(boundary::PeriodicExtended) = boundary.boundary
get_boundary_width(boundary::Periodic) = 0
function run_peristimulus_tests(raster_size, 
        boundary::AbstractBoundaryCondition, 
        lag_extents, 
        n_signals, signal_jitters,  noise_rate, 
        n_trials, n_tests, 
        postproc!, α
    )
    n_size, t_size = raster_size
    n_max_jitter, t_max_jitter = signal_jitters
    boundary_width = get_boundary_width(boundary)

    n0_range = (n_max_jitter+1):(n_size-n_max_jitter)
    t0_range = ((t_size ÷ 2) + t_max_jitter+boundary_width):(t_size - t_max_jitter-boundary_width)

    motif_class_range = 1:14

    test_results = DataFrame(mapreduce(vcat, motif_class_range, init=[]) do signal_motif_class_num
        l_detected_motifs = ThreadsX.map(1:n_tests) do test
            trials_epoch_tricorrs = jittered_trials_epochs(raster_size, boundary, lag_extents, signal_motif_class_num, n_signals, (n0_range, t0_range), signal_jitters, noise_rate, n_trials, postproc!)
            motif_significances = test_epoch_difference(trials_epoch_tricorrs).significance
            motif_significances .< α
        end
        proportion_detected_by_motif = mean(l_detected_motifs)
        ThreadsX.map(1:14) do detected_motif_num
            proportion_detected = proportion_detected_by_motif[detected_motif_num]
            (
                proportion_detected=proportion_detected,
                detected_motif=detected_motif_num,
                signal_motif=signal_motif_class_num,
                n_signals = n_signals,
                n_tests = n_tests,
                n_trials = n_trials
            )
        end
    end)
    return test_results
end


function test_peristimulus_difference(l_timeseries::AbstractVector{<:AbstractVector}, start, stop)
    peristimulus = vcat([timeseries[start:stop] for timeseries ∈ l_timeseries]...)
    nonstimulus = vcat([[timeseries[begin:start-1]; timeseries[stop+1:end]] for timeseries ∈ l_timeseries]...)
    effect_size = mean(peristimulus) - mean(nonstimulus)
    significance = pvalue(UnequalVarianceTTest(peristimulus, nonstimulus))
    return (effect_size=effect_size, significance=significance)
end

function safe_pvalue(a, b)
    if a == b
        return 1.
    else
        return pvalue(UnequalVarianceTTest(a, b))
    end
end

function test_epoch_difference(l_trials::AbstractVector{<:AbstractMatrix})
    # matrix of contributions
    n_motifs = size(first(l_trials),1)
    @assert n_motifs == 14
    @assert size(first(l_trials),2) == 2
    hypothesis_pairs = map(1:n_motifs) do i_motif
        stimulus = [trial[i_motif,2] for trial ∈ l_trials]
        nonstimulus = [trial[i_motif,1] for trial in l_trials]
        (stimulus, nonstimulus)
    end
    effect_size = hypothesis_pairs .|> (x) -> mean(x[1]) - mean(x[2])
    significance = [safe_pvalue(pair...) for pair in hypothesis_pairs]
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

function embedded_rand_motif(motif_class, n_size, t_size, n0_range::AbstractArray, t0_range::AbstractArray, n_max_jitter, t_max_jitter)
    raster = zeros(Bool, n_size, t_size)
    motif_coords = TripleCorrelations.rand_motif(motif_class, n0_range, t0_range, n_max_jitter, t_max_jitter)
    for motif_coord in motif_coords
        raster[CartesianIndex(motif_coord)] = 1
    end
    return raster
end

function embedded_rand_motif(motif_class, raster_size, base_node_ranges, signal_jitters)
    raster = zeros(Bool, raster_size...)
    motif_coords = TripleCorrelations.rand_motif(motif_class, base_node_ranges, signal_jitters)
    for motif_coord in motif_coords
        raster[CartesianIndex(motif_coord)] = 1
    end
    return raster
end

function calculate_trial_epochs(raster, boundary, lag_extents, epochs, postproc!)
    mapreduce(hcat, epochs) do epoch
        contributions = sequence_class_tricorr(raster[:,epoch], boundary, lag_extents)
        postproc!(contributions, contributions, raster[:,epoch])
        contributions
    end
end

@inline function view_slice_last(arr::AbstractArray{T,N}, dx) where {T,N}
    view(arr, ntuple(_ -> Colon(), N - 1)..., dx)
end

function fixed_noise_raster(raster, noise_rate, boundary::PeriodicExtended)
    dims = size(raster)
    n_bdry_ones = round(Int, noise_rate * prod(dims[1:end-1]) * boundary.boundary)
    n_meat_ones = round(Int, noise_rate * prod(dims[1:end-1]) * (dims[end] - 2(boundary.boundary)))
    n_meat_ones -= count(raster)
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

function fixed_noise_raster(raster, noise_rate, boundary)
    dims = size(raster)
    n_ones = round(Int, prod(dims) * noise_rate) - count(raster)
    noise_raster = zeros(Bool, dims...)
    noise_raster[1:n_ones] .= 1
    shuffle!(noise_raster)
    return noise_raster
end

function jittered_trials_epochs(raster_size::Tuple, 
        boundary::AbstractBoundaryCondition,
        lag_extents, 
        motif_class_num::Int, n_signals::Int, base_node_ranges, signal_jitters,
        noise_rate, n_trials,
        postproc!
    )
    n0_range, t0_range = base_node_ranges
    n_max_jitter, t_max_jitter = signal_jitters
    @assert t0_range[begin] > 1+t_max_jitter
    epochs = [1:t0_range[begin]-t_max_jitter,(t0_range[begin]-t_max_jitter+1):raster_size[end]] # HELP: is raster_size right?
    motif_class = offset_motif_numeral(motif_class_num)
    trials_epoch_tricorrs = map(1:n_trials) do _
        raster = embedded_rand_motif(motif_class, raster_size, base_node_ranges, signal_jitters)
        for _ in 1:(n_signals-1)
            raster .|= embedded_rand_motif(motif_class, raster_size, base_node_ranges, signal_jitters)
        end  
        for epoch in epochs
            epoch_raster = view_slice_last(raster, epoch)
            epoch_noise_raster = fixed_noise_raster(epoch_raster, noise_rate, boundary)
            epoch_raster .|= epoch_noise_raster
        end
        epoch_tricorrs = calculate_trial_epochs(raster, boundary, lag_extents, epochs, postproc!)
        epoch_tricorrs
    end
    return trials_epoch_tricorrs
end

