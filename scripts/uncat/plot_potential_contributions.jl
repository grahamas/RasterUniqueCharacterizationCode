using DrWatson; quickactivate(@__DIR__)

include(srcdir("main.jl"))

potential_fig = let N_neurons = 150, N_times = 150, max_lag = 7, transform=log10;

neuron_max_lag = time_max_lag = max_lag
raster = BitMatrix(ones(N_neurons, N_times))

fig = @time with_theme(bar_theme) do
    fig = Figure()
    fig[1,1] = potential_ax = CairoMakie.Axis(fig, yscale=transform)#, ygridvisible=true, yminorgridvisible=true, ygridcolor=:grey, yminorgridcolor=:lightgrey, yminorticks = IntervalsBetween(8))
    ret = plot_network_class_contributions!(potential_ax, raster, neuron_max_lag, time_max_lag)
    potential_ax.title[] = "Potential Contributions\n(max_lag = $(max_lag))"
    fig
end

save(plotsdir("potential_contributions_$(transform != identity ? string(transform) * "_" : "")lag$(max_lag)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")).png"), fig)

fig

end

