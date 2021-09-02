include(srcdir("main.jl"))
using Parameters


@inline tuple_join(x) = x
@inline tuple_join(x,y::AbstractString) = (x..., y)
@inline tuple_join(x::AbstractString,y) = (x, y...)
@inline tuple_join(x::Tuple,y::Tuple) = (x..., y...)
@inline tuple_join(x,y,z...) = (x..., tuple_join(y,z...)...)

function calc_counts_coord(n1, t1, n2, t2, accumulated_lag_contrib, neuron_max_lag, time_max_lag, raster_offset=1)
    n1 += neuron_max_lag + 1; n2 += neuron_max_lag + 1
    t1 += time_max_lag + 1; t2 += time_max_lag + 1

    major_grid_n = (neuron_max_lag * 2 + 1) * (n1 - 1) + 1
    major_grid_t = (time_max_lag * 2 + 1) * (t1 - 1) + 1
    
    grid_n = major_grid_n + n2 - 1
    grid_t = major_grid_t + t2 - 1
    (grid_t, accumulated_lag_contrib, grid_n+raster_offset)
end


function calc_raster_coord(n0, t0, nl, tl, (n_raster_neurons, n_raster_times), neuron_max_lag, time_max_lag, max_contrib, raster_depth=0)
    n_neuron_lags, n_time_lags = (neuron_max_lag, time_max_lag) .* 2 .+ 1
    time_lag_bins = (n_time_lags)^2-1
    time_multiple = (time_lag_bins) / (n_raster_times-1) 
    neuron_multiple = max_contrib / (n_raster_neurons-1) # If neuron=1, want coord=0; if neuron=9, want coord=50
    ((t0+tl-1)*time_multiple, (n0+nl-1)*neuron_multiple, raster_depth)
end

function calc_accumulating_row(n0, t0, n1, t1, n2, t2, i_step, n_steps, raster, accumulated_lag_contrib, neuron_max_lag, time_max_lag, max_contrib=10)
    raster_coord0 = calc_raster_coord(n0, t0, 0, 0, size(raster), neuron_max_lag, time_max_lag, max_contrib)
    contrib = raster[n0,t0] * raster[n0+n1,t0+t1] * raster[n0+n2,t0+t2]
    if n1 == t1 == n2 == t2 == 0
        @show n0 t0 contrib
    end
    if contrib != 0
        i_raster_appearance = 2(i_step-1) + 1
        i_counts_appearance = i_raster_appearance + 1
        motif_num = info_flow_classify_lag_motif_class(n1, n2, t1, t2)
        motif_numeral = roman_encode(motif_num)
        accumulated_contrib = accumulated_lag_contrib[n1, t1, n2, t2] += contrib
        counts_coord = calc_counts_coord(n1, t1, n2, t2, accumulated_contrib, neuron_max_lag, time_max_lag)
        raster_coords = [raster_coord0, 
            calc_raster_coord(n0, t0, n1, t1, size(raster), neuron_max_lag, time_max_lag, max_contrib), calc_raster_coord(n0, t0, n2, t2, size(raster), neuron_max_lag, time_max_lag, max_contrib)]
        error_coord = (-1, -1, -1)
        node_rows = map(1:3) do i_node
            tuple_join(n0, t0, n1, t1, n2, t2, i_node-1, [error_coord for _ ∈ 1:(i_raster_appearance-1)]..., raster_coords[i_node], [counts_coord for _ ∈ i_counts_appearance:2n_steps]..., i_step, motif_numeral)
        end
        return node_rows
    end
    return []
end


let fn = datadir("flow_illustration","manual_raster_7x7.csv"), neuron_max_lag=2,
    time_max_lag=2
    raster = load_raster(fn)

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

