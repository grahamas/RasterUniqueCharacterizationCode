function triple_correlation_class_contributions!(class_contribution::Vector, raster::BitMatrix, neuron_max_lag, time_max_lag, lags_classifier::Function)
    neuron_lag_range = -(neuron_max_lag):(neuron_max_lag)        
    time_lag_range = -(time_max_lag):(time_max_lag)

    (N_neurons, N_times) = size(raster)
    time_range = (1-minimum(time_lag_range)):(N_times-maximum(time_lag_range))
    neuron_range = (1-minimum(neuron_lag_range)):(N_neurons-maximum(neuron_lag_range))

    class_contribution .= 0

    @turbo for n1 ∈ neuron_lag_range, n2 ∈ neuron_lag_range, 
            t1 ∈ time_lag_range, t2 ∈ time_lag_range
        class = lags_classifier(n1, n2, t1, t2)
        contribution = 0
        for i_neuron ∈ neuron_range, i_time ∈ time_range #(i_neuron, i_time) ∈ IterTools.product(neuron_range, time_range)
            contribution += raster[i_neuron, i_time] * raster[i_neuron+n1,i_time+t1] * raster[i_neuron+n2,i_time+t2]
        end
        class_contribution[class] += contribution
    end
    return class_contribution ./ calculate_scaling_factor(raster, (neuron_max_lag, time_max_lag))
end

function triple_correlation_network_classifications(raster::BitMatrix, neuron_max_lag, time_max_lag)
    N_network_classifications = 14
    network_class_contributions = Array{Int}(undef, N_network_classifications)
    lags_classifier = info_flow_classify_lag_motif_class

    triple_correlation_class_contributions!(network_class_contributions, raster, neuron_max_lag, time_max_lag, lags_classifier)
end

function calculate_unscaled_triple_correlation!(correlation::OffsetArray{T,4}, raster::BitMatrix, λ_max) where T
    calculate_unscaled_triple_correlation!(parent(correlation), raster, λ_max)
end

function calculate_unscaled_triple_correlation!(correlation::Array{T,4}, raster::BitMatrix, λ_max) where T
    single_λ_ranges = [-l:l for l ∈ λ_max]
    neuron_lag_range, time_lag_range = single_λ_ranges


    (N_neurons, N_times) = size(raster)
    time_range = (1-minimum(time_lag_range)):(N_times-maximum(time_lag_range))
    neuron_range = (1-minimum(neuron_lag_range)):(N_neurons-maximum(neuron_lag_range))

    @turbo for i_n1 ∈ eachindex(neuron_lag_range), i_n2 ∈ eachindex(neuron_lag_range),  i_t1 ∈ eachindex(time_lag_range), i_t2 ∈ eachindex(time_lag_range)
        n1, n2, t1, t2 = neuron_lag_range[i_n1], neuron_lag_range[i_n2], 
time_lag_range[i_t1], time_lag_range[i_t2]
        accum = 0
        for i_neuron ∈ neuron_range, i_time ∈ time_range
            accum += raster[i_neuron, i_time] * raster[i_neuron+n1,i_time+t1] * raster[i_neuron+n2,i_time+t2]
        end
        correlation[i_n1,i_n2,i_t1,i_t2] = accum
    end
end

function calculate_unscaled_triple_correlation(raster::BitMatrix, λ_max::Tuple)
    λ_ranges = Tuple(repeat([-l:l for l ∈ λ_max], outer=(2,)))
    correlation = zeros(Float64, λ_ranges)

    calculate_unscaled_triple_correlation!(correlation, raster, λ_max)
    return correlation
end

function calculate_scaled_triple_correlation(raster::BitMatrix, λ_max::Tuple)
    unscaled_correlation = calculate_unscaled_triple_correlation(raster, λ_max)
    unscaled_correlation ./= calculate_scaling_factor(raster, λ_max)
end

function calculate_scaling_factor(raster, λ_max)
    N, T = size(raster)
    (T - λ_max[1] + 1) * (N - λ_max[2] + 1)
end

_hi_bound(l1, l2, M) = M - max(0, l1, l2)
_lo_bound(l1, l2) = 1 + min(0, l1, l2)

function n_triplets(n1, t1, n2, t2, N, T)
    (_hi_bound(n1, n2, N) - _lo_bound(n1, n2) + 1) *
    (_hi_bound(t1, t2, T) - _lo_bound(t1, t2) + 1)
end