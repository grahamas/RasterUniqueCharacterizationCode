
function load_raster(fn)
    return CSV.read(fn, BitArray ∘ Tables.matrix, delim=",", header=false)
end

function write_raster(fn, arr)
    tbl = Tables.table(arr)
    CSV.write(fn, tbl, header=false)
end

function generate_random_raster(raster_size::Tuple, spike_prob=0.1)
    putative_inputs = rand(Float64, raster_size)
    spikes = BitArray(putative_inputs .>= spike_prob)
end

