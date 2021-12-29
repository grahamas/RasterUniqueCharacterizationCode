function make_neuron_raster(dims...)
    arr = BitArray(undef, dims...)
    arr .= 0
    make_neuron_raster(arr)
end
function make_neuron_raster(arr::BitMatrix)
    dims = size(arr)
    OffsetArray(arr, [1:dim for dim in dims[begin:end-1]]..., 0:(dims[end]-1))
end

function load_raster(fn)
    arr = CSV.read(fn, BitArray ∘ Tables.matrix, delim=",", header=false)
    return OffsetArray(arr, 1:size(arr,1), 0:(size(arr,2)-1))
end

function write_raster(fn, arr)
    tbl = Tables.table(arr)
    CSV.write(fn, tbl, header=false)
end

function write_raster_coordinates(fn, arr::OffsetArray)
    coords = Tuple.(findall(arr))
    time_zeroed_coords = [(x[1], x[2]) for x in coords]
    df = DataFrame(time_zeroed_coords)
    rename!(df, ["Neuron", "Time"])
    CSV.write(fn, df)
end
function write_raster_coordinates(fn, arr::BitMatrix)
    write_raster_coordinates(fn, make_neuron_raster(arr))
end

function generate_random_raster(raster_size::Tuple, spike_prob=0.1)
    putative_inputs = rand(Float64, raster_size)
    BitArray(putative_inputs .<= spike_prob)
end

function write_spike_trains(fn, arr)
    df = DataFrame(Neuron=Int[], Time=Int[], Spike=Int[])
    for i ∈ CartesianIndices(arr)
        spike = arr[i]
        i_neuron, i_time = Tuple(i)
        push!(df, (Neuron=i_neuron, Time=i_time, Spike=Int(spike)))
    end
    CSV.write(fn, df)
end

function base_nodes_raster(raster::OffsetArray, neuron_max_lag, time_max_lag)
    neuron_lag_range = -(neuron_max_lag):(neuron_max_lag)        
    time_lag_range = -(time_max_lag):(time_max_lag)

    (N_neurons, N_times) = size(raster)
    time_range = (0-minimum(time_lag_range)):(N_times-1-maximum(time_lag_range))
    neuron_range = (1-minimum(neuron_lag_range)):(N_neurons-maximum(neuron_lag_range))

    base_nodes_raster = copy(raster)
    base_nodes_raster .= zero(eltype(raster))
    for i_neuron ∈  neuron_range, i_time ∈ time_range
        base_nodes_raster[i_neuron, i_time] = one(eltype(raster))
    end
    return base_nodes_raster
end