@quickactivate "RasterUniqueCharacterizationCode"

using TripleCorrelations
using UnicodePlots

raster_size = (1500, 2000)
lag_extents = (5,7)

raster = zeros(Bool, raster_size...)

expectations = map(100:100:5000000) do n_spikes
    raster[1:n_spikes] .= 1
    expectation_conditioned_on_spike_count(raster, Periodic(), lag_extents)
end

lineplot((getindex.(expectations, 2) .* getindex.(expectations, 1)), getindex.(expectations, 3))