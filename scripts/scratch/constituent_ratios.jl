@quickactivate "RasterUniqueCharacterizationCode"

using TripleCorrelations
using UnicodePlots

raster_size = (1500, 2000)
lag_extents = (5,7)

expectations = map(100:100:500000) do n_spikes
    expectation_conditioned_on_spike_count(n_spikes, raster_size, lag_extents)
end

lineplot((getindex.(expectations, 4) .* getindex.(expectations, 1)), getindex.(expectations, 5))