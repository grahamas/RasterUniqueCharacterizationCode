

include(srcdir("main.jl"))

const InfoFlowRow = @NamedTuple{spike_motif::String, n1_num::Int, n2_num::Int, t1_num::Int, t2_num::Int, motif_class_num::Int}
lag_1(ifr::InfoFlowRow) = TernarizedLag((ifr.n1_num, ifr.t1_num))
lag_2(ifr::InfoFlowRow) = TernarizedLag((ifr.n2_num, ifr.t2_num))
#lag_triplet(ifr::InfoFlowRow) = LagTriplet(lag_1(ifr), lag_2(ifr))
lag_motif(ifr::InfoFlowRow) = LagMotif(lag_1(ifr), lag_2(ifr))


let raster_size = (150, 150),
        λ_max = (7, 7),
        freq = 0.12, θ=0, noise_amplitude = 0, signal_amplitude=1,
        spike_rate=0.5;
    lag_classes = generate_all_lag_motif_classes()
    info_flow_classes = info_flow_classify_lag_motif_classes(lag_classes)

    random_raster = generate_random_raster(raster_size, spike_rate)
    potential_raster = BitArray(ones(Bool, raster_size))
    raster_coords = IterTools.product(1:raster_size[1], 1:raster_size[2])
    raster_signal = map(raster_coords) do (i_neuron, i_time)
        (sin(2π*freq*i_time + θ) + 1) / 2
    end
    raster_noise = noise_amplitude * rand(Float64, raster_size)
    structured_raster = (raster_signal + raster_noise) .> spike_rate
    
    random_triple_correlation = triple_correlation(random_raster, λ_max)
    potential_triple_correlation = triple_correlation(potential_raster, λ_max)
    structured_triple_correlation = triple_correlation(structured_raster, λ_max)

    random_lag_motifs = reduce_to_lag_motifs(random_triple_correlation)
    potential_lag_motifs = reduce_to_lag_motifs(potential_triple_correlation)
    structured_lag_motifs = reduce_to_lag_motifs(structured_triple_correlation)


    df = map(info_flow_classes) do row
        named_row = name_info_flow_row(row)

        random=random_lag_motifs[lag_motif(named_row)]
        potential=potential_lag_motifs[lag_motif(named_row)]
        structured=structured_lag_motifs[lag_motif(named_row)]

        merge(named_row, 
            (
                random=random,
                potential=potential,
                structured=structured,
                random_ratio=random./potential,
                structured_pratio=structured./potential,
                structured_rratio=structured./random
             )
        )
    end |> DataFrame
    CSV.write(datadir("flow", "prevalences.csv"), df)
    write_raster_coordinates(datadir("flow", "random_raster_coords.csv"), random_raster)
    write_raster_coordinates(datadir("flow", "potential_raster_coords.csv"), potential_raster)
    write_raster_coordinates(datadir("flow", "structured_raster_coords.csv"), structured_raster)
end
