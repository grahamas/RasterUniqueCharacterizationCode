@quickactivate "RasterUniqueCharacterizationCode"

using TripleCorrelations, Random, Statistics

#defines motif_examples dict
include("../../TripleCorrelations/test/data/motif_examples.jl")
include("../../TripleCorrelations/test/src/helpers.jl")

function an_ratio(raster, n_lag, t_lag; noise_calcs=1)
    raw_contribution = sequence_class_tricorr(raster, n_lag, t_lag)
    noise_contribution = sum(sequence_class_tricorr(shuffle(raster, n_lag, t_lag)) for _ in 1:noise_calcs) ./ noise_calcs
    raw_contribution ./ noise_contribution
end

function an_timeseries(raster, n_lag, t_lag, step; t_window=2t_lag+1)
    N,T = size(raster)
    map(1:step:(T-t_window)) do t_start
        an_ratio(raster[:,t_start:t_start+t_window], n_lag, t_lag)
    end
end


an_contributions = let n_pad = 5, t_pad = 30, 
    n_reps = 1, t_reps = 1,
    n_lag = 6, t_lag = 10,
    snippets=10, step=2,
    noise_rate = 0.1;

signal_raster = repeat_padded_motif("III", n_pad, t_pad, n_reps, t_reps)
noise_raster = rand(size(signal_raster)...) .< noise_rate

raster = max.(signal_raster .+ noise_raster, 1)

total_contributions = sequence_class_tricorr(raster, n_lag, t_lag)
noise_contributions = sequence_class_tricorr(shuffle(raster), n_lag, t_lag)
an_contributions = total_contributions ./ noise_contributions
an_contributions
end

