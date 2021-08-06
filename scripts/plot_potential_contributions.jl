using DrWatson; quickactivate(@__DIR__)

include(srcdir("main.jl"))

potential_fig = let N_neurons = 150, N_times = 150, window_size = 15;

neuron_window_size = time_window_size = 15
raster = BitMatrix(ones(N_neurons, N_times))

with_theme(bar_theme) do
    fig = Figure()
    fig[1,1] = potential_ax = CairoMakie.Axis(fig)
    ret = plot_network_class_contributions!(potential_ax, raster, neuron_window_size, time_window_size)
    potential_ax.title[] = "Potential Contributions"
    fig
end

end



save(plotsdir("potential_contributions.png"), potential_fig)