using DrWatson; quickactivate(@__DIR__)

include(srcdir("main.jl"))

potential_fig = let N_neurons = 150, N_times = 150, window_size = 15, transform=log10;

neuron_window_size = time_window_size = window_size
raster = BitMatrix(ones(N_neurons, N_times))

fig = with_theme(bar_theme) do
    fig = Figure()
    fig[1,1] = potential_ax = CairoMakie.Axis(fig, yscale=transform)#, ygridvisible=true, yminorgridvisible=true, ygridcolor=:grey, yminorgridcolor=:lightgrey, yminorticks = IntervalsBetween(8))
    fillto = if transform == log10
        potential_ax.ytickformat[] = Makie.automatic
        1.
    else
        0.
    end
    ret = plot_network_class_contributions!(potential_ax, raster, neuron_window_size, time_window_size; fillto=fillto)
    potential_ax.title[] = "Potential Contributions\n(window size = $(window_size))"
    fig
end

save(plotsdir("potential_contributions_$(transform != identity ? string(transform) * "_" : "")lag$(window_size)_$(Dates.now()).png"), fig)

fig

end

