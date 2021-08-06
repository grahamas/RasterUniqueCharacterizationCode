
function load_raster(fn)
    return CSV.read(fn, BitArray ∘ Tables.matrix, delim=",", header=false)
end

function write_raster(fn, arr)
    tbl = Tables.table(arr)
    CSV.write(fn, tbl, header=false)
end

function random_raster(N_neurons, N_times, spike_prob=0.1)
    putative_inputs = rand(N_neurons, N_times)
    spikes = BitArray(putative_inputs .>= spike_prob)
end