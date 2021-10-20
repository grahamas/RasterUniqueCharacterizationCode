
let neuron_max_lag=2,
        time_max_lag=2,
        N_neurons = 30, N_times = 30,
        firing_threshold = 0.7;

    all_ones_raster = make_neuron_raster(N_neurons, N_times)
    all_ones_raster .= 1
    all_ones_contributions = triple_correlation_network_classifications(all_ones_raster, neuron_max_lag, time_max_lag)

    random_raster = make_neuron_raster(N_neurons, N_times)
    for i ∈ eachindex(random_raster)
        random_raster[i] = rand() >= firing_threshold
    end
    random_contributions = triple_correlation_network_classifications(random_raster, neuron_max_lag, time_max_lag)

    banded_contributions, banded_raster = let freq = 0.12, θ=0, noise_amplitude = 0, signal_amplitude=1;
        raster_coords = IterTools.product(1:N_neurons, 1:N_times)
        raster_signal = map(raster_coords) do (i_neuron, i_time)
            (sin(2π*freq*i_time + θ) + 1) / 2
        end
        raster_noise = noise_amplitude * rand(N_neurons, N_times)
        raster = (raster_signal + raster_noise) .> firing_threshold

        contributions = triple_correlation_network_classifications(raster, neuron_max_lag, time_max_lag)
        (contributions, raster)
    end

    shifting_contributions, shifting_raster = let freq = 0.12, θ=-2π, noise_amplitude = 0;
        raster_coords = IterTools.product(1:N_neurons, 1:N_times)
        raster_signal = map(raster_coords) do (i_neuron, i_time)
            (sin(2π*freq*i_time + θ*i_neuron/N_neurons) + 1) / 2
        end
        raster_noise = noise_amplitude * rand(N_neurons, N_times)
        raster = (raster_signal + raster_noise) .> firing_threshold
    
        (triple_correlation_network_classifications(raster, neuron_max_lag, time_max_lag), raster)
    end

    tonic_contributions, tonic_raster = let freq = 0.12, θ=-0, noise_amplitude = 0;
        raster_coords = IterTools.product(1:N_neurons, 1:N_times)
        raster_signal = map(raster_coords) do (i_neuron, i_time)
            (sin(2π*freq*i_time + θ*i_neuron/N_neurons) + 1) / 2 + (i_neuron % 10 == 0 ? 1 : 0)
        end
        raster_noise = noise_amplitude * rand(N_neurons, N_times)
        raster = (raster_signal + raster_noise) .> firing_threshold
    
        (triple_correlation_network_classifications(raster, neuron_max_lag, time_max_lag), raster)
    end

    low_noise_banded_contributions, low_noise_banded_raster = let freq = 0.12, θ=0, noise_amplitude=0.5, signal_amplitude=0.5;
        raster_coords = IterTools.product(1:N_neurons, 1:N_times)
        raster_signal = map(raster_coords) do (i_neuron, i_time)
            (sin(2π*freq*i_time + θ) + 1) / 2
        end
        raster_noise = rand(N_neurons, N_times)
        raster = (signal_amplitude .* raster_signal + noise_amplitude .* raster_noise) .> firing_threshold
    
        (triple_correlation_network_classifications(raster, neuron_max_lag, time_max_lag), raster)
    end

    mid_noise_banded_contributions, mid_noise_banded_raster = let freq = 0.12, θ=0, noise_amplitude=2//3, signal_amplitude=1//3;
        raster_coords = IterTools.product(1:N_neurons, 1:N_times)
        raster_signal = map(raster_coords) do (i_neuron, i_time)
            (sin(2π*freq*i_time + θ) + 1) / 2
        end
        raster_noise = rand(N_neurons, N_times)
        raster = (signal_amplitude .* raster_signal + noise_amplitude .* raster_noise) .> firing_threshold
    
        (triple_correlation_network_classifications(raster, neuron_max_lag, time_max_lag), raster)
    end
    
    high_noise_banded_contributions, high_noise_banded_raster= let freq = 0.04, θ=0, noise_amplitude=1;
        raster_coords = IterTools.product(1:N_neurons, 1:N_times)
        raster_signal = map(raster_coords) do (i_neuron, i_time)
            (sin(2π*freq*i_time + θ) + 1) / 2
        end
        raster_noise = noise_amplitude * rand(N_neurons, N_times)
        raster = (raster_noise) .> firing_threshold
    
        (triple_correlation_network_classifications(raster, neuron_max_lag, time_max_lag), raster)
    end

    all_ones_tc = triple_correlation(all_ones_raster, (neuron_max_lag, time_max_lag))

    df = DataFrame()
    df.Motifs = roman_encode.(1:14) |> collect
    
    df.all_ones = all_ones_contributions
    write_raster_coordinates(datadir("flow_illustration", "all_ones_raster_30x30.csv"), all_ones_raster)

    df.random = random_contributions
    df.random_of_possible = random_contributions ./ all_ones_contributions
    write_raster_coordinates(datadir("flow_illustration", "random_raster_30x30.csv"), random_raster)

    df.banded = banded_contributions
    df.banded_of_possible = banded_contributions ./ all_ones_contributions
    df.banded_of_random = banded_contributions ./ random_contributions
    write_raster_coordinates(datadir("flow_illustration", "banded_raster_30x30.csv"), banded_raster)

    df.shifting = shifting_contributions
    df.shifting_of_possible = shifting_contributions ./ all_ones_contributions
    df.shifting_of_random = shifting_contributions ./ random_contributions
    write_raster_coordinates(datadir("flow_illustration", "shifting_raster_30x30.csv"), shifting_raster)

    df.tonic = tonic_contributions
    df.tonic_of_possible = tonic_contributions ./ all_ones_contributions
    df.tonic_of_random = tonic_contributions ./ random_contributions
    write_raster_coordinates(datadir("flow_illustration", "tonic_raster_30x30.csv"), tonic_raster)

    df.low_noise_banded = low_noise_banded_contributions
    df.low_noise_banded_of_possible = low_noise_banded_contributions ./ all_ones_contributions
    df.low_noise_banded_of_random = low_noise_banded_contributions ./ random_contributions
    write_raster_coordinates(datadir("flow_illustration", "low_noise_banded_raster_30x30.csv"), low_noise_banded_raster)

    df.mid_noise_banded = mid_noise_banded_contributions
    df.mid_noise_banded_of_possible = mid_noise_banded_contributions ./ all_ones_contributions
    df.mid_noise_banded_of_random = mid_noise_banded_contributions ./ random_contributions
    write_raster_coordinates(datadir("flow_illustration", "mid_noise_banded_raster_30x30.csv"), mid_noise_banded_raster)

    df.high_noise_banded = high_noise_banded_contributions
    df.high_noise_banded_of_possible = high_noise_banded_contributions ./ all_ones_contributions
    df.high_noise_banded_of_random = high_noise_banded_contributions ./ random_contributions
    write_raster_coordinates(datadir("flow_illustration", "high_noise_banded_raster_30x30.csv"), high_noise_banded_raster)

    CSV.write(datadir("flow_illustration", "results_30x30.csv"), df)

    # nts = map(rows) do row
    #     NamedTuple{header_tuple}(row)
    # end
    # df = DataFrame(nts)

    # base_nodes_arr = base_nodes_raster(raster, neuron_max_lag, time_max_lag)
    # write_raster_coordinates(datadir("flow_illustration", "base_nodes_raster_7x7.csv"), base_nodes_arr)

    # write_spike_trains(datadir("flow_illustration", "spike_trains_7x7.csv"), raster)
    # write_raster_coordinates(datadir("flow_illustration", "raster_7x7.csv"), raster)
    # CSV.write(datadir("flow_illustration", "illustrated_7x7.csv"), df)
end

