let neuron_max_lag=2,
    time_max_lag=2
    raster = load_raster(datadir("flow_illustration","manual_raster_7x7.csv"))

    neuron_lag_range = -(neuron_max_lag):(neuron_max_lag)        
    time_lag_range = -(time_max_lag):(time_max_lag)

    (N_neurons, N_times) = size(raster)
    time_range = (0-minimum(time_lag_range)):(N_times-1-maximum(time_lag_range))
    neuron_range = (1-minimum(neuron_lag_range)):(N_neurons-maximum(neuron_lag_range))

    accumulated_lag_contrib = zeros(Int, neuron_lag_range, time_lag_range, neuron_lag_range, time_lag_range)

    @warn "Looping $(length(neuron_range)*length(time_range)*length(neuron_lag_range)^2*length(time_lag_range)^2) times..."
    
    step_points = [
        (location=(3,2), lag=(2,0,1,1)),
        (location=(3,2), lag=(0,0,0,0)),
        (location=(5,2), lag=(-1,1,-2,0))
    ]
    post_step_points = length(step_points) + 1
    header = vcat(
        ["n", "t", "n1", "t1", "n2", "t2", "node"],
        [
            vcat(
                [
                    ["$(name)_x", "$(name)_y", "$(name)_z"] 
                    for name ∈ ["raster_$(SP)", "counts_$(SP)"]
                ]...
            ) for SP in 1:post_step_points
        ]...
    )
    push!(header, "appearance")
    push!(header, "Motif")
    header_tuple = Symbol.(tuple(header...))
    rows = NTuple{length(header_tuple),Union{Float64,String}}[]
    for (step_num, step_point) ∈ enumerate(step_points)
        @unpack location, lag = step_point
        n1, t1, n2, t2 = lag
        raster[location...]*raster[(location .+ (n1,t1))...]*raster[(location .+ (n2,t2))...] == 0 && @warn "Non-extant triplet: $step_point"
        triplet_rows = calc_accumulating_row(location..., n1, t1, n2, t2,
                            step_num, post_step_points, raster, 
                            accumulated_lag_contrib,
                            neuron_max_lag, time_max_lag)
        rows = vcat(rows, triplet_rows)
    end
    @show neuron_range, time_range
    for i_neuron ∈ neuron_range, i_time ∈ time_range
        # sub-optimal loop nesting for storytelling purposes
        if raster[i_neuron, i_time] != 0
            row_triples = [calc_accumulating_row(i_neuron, i_time, n1, t1, n2, t2, 
                                        post_step_points, post_step_points,
                                        raster, accumulated_lag_contrib,
                                        neuron_max_lag, time_max_lag
                                    ) 
                            for n1 ∈ neuron_lag_range, n2 ∈ neuron_lag_range,
                                t1 ∈ time_lag_range, t2 ∈ time_lag_range
                            if (location=(i_neuron, i_time), lag=(n1, t1, n2, t2)) ∉ step_points
            ]
            rows = vcat(rows, row_triples...)
        end
    end

    nts = map(rows) do row
        NamedTuple{header_tuple}(row)
    end
    df = DataFrame(nts)

    base_nodes_arr = base_nodes_raster(raster, neuron_max_lag, time_max_lag)
    write_raster_coordinates(datadir("flow_illustration", "base_nodes_raster_7x7.csv"), base_nodes_arr)

    write_spike_trains(datadir("flow_illustration", "spike_trains_7x7.csv"), raster)
    write_raster_coordinates(datadir("flow_illustration", "raster_7x7.csv"), raster)
    CSV.write(datadir("flow_illustration", "illustrated_7x7.csv"), df)
end

