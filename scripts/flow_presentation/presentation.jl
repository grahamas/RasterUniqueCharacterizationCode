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
    raster_neuron_lim = 8
    raster_time_lim = 8
    n_neuron_lags, n_time_lags = (neuron_max_lag, time_max_lag) .* 2 .+ 1
    time_lag_bin_max = (n_time_lags)^2
    time_multiple = (time_lag_bin_max) / (raster_time_lim) 
    neuron_multiple = max_contrib / (raster_neuron_lim) # If neuron=1, want coord=0; if neuron=9, want coord=50
    ((t0+tl)*time_multiple, (n0+nl)*neuron_multiple, raster_depth)
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


include("illustration_7x7.jl")
include("results_30x30.jl")