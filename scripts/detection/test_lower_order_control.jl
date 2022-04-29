@quickactivate "RasterUniqueCharacterizationCode"

using TripleCorrelations
using HypothesisTests, Random, Statistics
using ProgressMeter
using CairoMakie; plot_ext = "png"
using Dates
using Pkg
using Base.Threads
using DataFrames, AlgebraOfGraphics, Colors
using JLD2

include(srcdir("roman_encode.jl"))
include(srcdir("peristimulus_testing.jl"))

(motif_I_2x, motif_II_x) = let x = 20,
    n_size = 32, t_size = 60,
    n_max_jitter = 3, t_max_jitter = 2,
    lag_extents = (12, 10);

n0_range = (n_max_jitter+1):(n_size-n_max_jitter)
t0_range = (1 + t_max_jitter):(t_size - t_max_jitter)

raster_I = embedded_rand_motif("I", n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter)
for i ∈ 2:2x
    raster_I .|= embedded_rand_motif("I", n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter)
end

raster_II = embedded_rand_motif("II", n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter)
for i ∈ 2:x
    raster_II .|= embedded_rand_motif("II", n_size, t_size, n0_range, t0_range, n_max_jitter, t_max_jitter)
end

contributions_I = sequence_class_tricorr(raster_I, Periodic(), lag_extents) .* TripleCorrelations.calculate_scaling_factor(raster_I, Periodic())
contributions_II = sequence_class_tricorr(raster_II, Periodic(), lag_extents) .* TripleCorrelations.calculate_scaling_factor(raster_II, Periodic())

rate_control_I = expectation_conditioned_on_spike_count(count(raster_I), size(raster_I), lag_extents)
rate_control_II = expectation_conditioned_on_spike_count(count(raster_II), size(raster_II), lag_extents)

lower_order_control_I = expectation_conditioned_on_lower_orders(contributions_I, count(raster_I), size(raster_I), lag_extents)
lower_order_control_II = expectation_conditioned_on_lower_orders(contributions_II, count(raster_II), size(raster_II), lag_extents)

((raster_I, contributions_I, contributions_I ./ rate_control_I, contributions_I ./ lower_order_control_I), (raster_II, contributions_II, contributions_II ./ rate_control_II, contributions_II ./ lower_order_control_II))

end