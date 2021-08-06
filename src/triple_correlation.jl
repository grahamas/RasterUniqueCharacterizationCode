function triple_correlation_class_contributions!(class_contribution::Vector, raster::BitMatrix, neuron_window, time_window, lags_classifier::Function)
    neuron_lag_range = -(neuron_window ÷ 2):(neuron_window ÷ 2)        
    time_lag_range = -(time_window ÷ 2):(time_window ÷ 2)

    (N_neurons, N_times) = size(raster)
    time_range = (1-minimum(time_lag_range)):(N_times-maximum(time_lag_range))
    neuron_range = (1-minimum(neuron_lag_range)):(N_neurons-maximum(neuron_lag_range))

    for n1 ∈ neuron_lag_range, n2 ∈ neuron_lag_range, 
            t1 ∈ time_lag_range, t2 ∈ time_lag_range
        class = lags_classifier(n1, n2, t1, t2)
        contribution = 0
        @turbo for i_neuron ∈ neuron_range, i_time ∈ time_range #(i_neuron, i_time) ∈ IterTools.product(neuron_range, time_range)
            contribution += raster[i_neuron, i_time] * raster[i_neuron+n1,i_time+t1] * raster[i_neuron+n2,i_time+t2]
        end
        class_contribution[class] += contribution
    end
    return class_contribution
end

function triple_correlation_network_classifications(raster::BitMatrix, neuron_window, time_window)
    N_network_classifications = 14
    network_class_contributions = zeros(Int, N_network_classifications)
    lags_classifier = classify_network_motif

    triple_correlation_class_contributions!(network_class_contributions, raster, neuron_window, time_window, lags_classifier)
end
