
include(srcdir("main.jl"))


TT_fig, NN_fig = let freq=0.1, θ=-2π, N_neurons=150, N_times=150,
    neuron_window=31, time_window=31;

raster_coords = IterTools.product(1:N_neurons, 1:N_times)
raster_signal = map(raster_coords) do (i_neuron, i_time)
    (sin(2π*freq*i_time + θ*i_neuron/N_neurons) + 1) / 2
end
raster = BitArray(raster_signal .> 0.5)
triple_correlation = calculate_unscaled_triple_correlation(raster, neuron_window, time_window)

TT = dropdims(mean(triple_correlation, dims=(1,2)), dims=(1,2))
NN = dropdims(mean(triple_correlation, dims=(3,4)), dims=(3,4))

TT_fig = Figure()
TT_fig[1,1] = TT_ax = CairoMakie.Axis3(TT_fig)
surface!(TT_ax, collect.(axes(TT))..., parent(TT))

NN_fig = Figure()
NN_fig[1,1] = NN_ax = CairoMakie.Axis3(NN_fig)
surface!(NN_ax, collect.(axes(NN))..., parent(NN))

(TT_fig, NN_fig)

end
