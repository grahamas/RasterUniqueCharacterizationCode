
include(srcdir("main.jl"))

fig3 = with_theme(bar_theme) do

fig3 = Figure(resolution=(1600, 2400))

let N_neurons = 150, N_times = 150,
        lag_extents = (14, 14),
        boundary = Periodic();
all_ones_raster = Array{Bool}(ones(N_neurons, N_times))
potential_contributions = sequence_class_tricorr_unrolled(
    all_ones_raster, boundary, lag_extents
)

subfig3A = let freq = 0.12, θ=0, noise_amplitude = 0, signal_amplitude=1;
    raster_coords = IterTools.product(1:N_neurons, 1:N_times)
    raster_signal = map(raster_coords) do (i_neuron, i_time)
        (sin(2π*freq*i_time + θ) + 1) / 2
    end
    raster_noise = noise_amplitude * rand(N_neurons, N_times)
    raster = (raster_signal + raster_noise) .> 0.5

   # contributions_axis = CairoMakie.Axis(fig3)
    relative_contributions_axis = CairoMakie.Axis(fig3)
    plot_relative_network_class_contributions!(relative_contributions_axis, raster, boundary, potential_contributions, lag_extents)
    raster_axis = CairoMakie.Axis(fig3)    
    # plot_network_class_contributions!(contributions_axis, raster, lag_extents)
    # hideydecorations!(contributions_axis, ticks=false, ticklabels=false)
    # hidexdecorations!(contributions_axis)
    hideydecorations!(relative_contributions_axis, ticks=false, ticklabels=false)
    hidexdecorations!(relative_contributions_axis)

    plot_raster!(raster_axis, raster)

    sublayout = GridLayout()
    title = Label(fig3, "freq = $freq AU", tellheight=true, tellwidth=false)
    sublayout[:v] = [title, raster_axis, relative_contributions_axis]

    sublayout
end

fig3[1,1] = subfig3A

fig3[1,2] = subfig3B = let freq = 0.12, θ=-2π, noise_amplitude = 0;
    raster_coords = IterTools.product(1:N_neurons, 1:N_times)
    raster_signal = map(raster_coords) do (i_neuron, i_time)
        (sin(2π*freq*i_time + θ*i_neuron/N_neurons) + 1) / 2
    end
    raster_noise = noise_amplitude * rand(N_neurons, N_times)
    raster = (raster_signal + raster_noise) .> 0.5

   # contributions_axis = CairoMakie.Axis(fig3)
    relative_contributions_axis = CairoMakie.Axis(fig3)
    plot_relative_network_class_contributions!(relative_contributions_axis, raster, boundary, potential_contributions, lag_extents)
    raster_axis = CairoMakie.Axis(fig3)
    sublayout = GridLayout()
    title = Label(fig3, "phase shifting", tellheight=true, tellwidth=false)
    sublayout[:v] = [title, raster_axis, relative_contributions_axis]
    # plot_network_class_contributions!(contributions_axis, raster, lag_extents)
    # hideydecorations!(contributions_axis, ticks=false, ticklabels=false)
    # hidexdecorations!(contributions_axis)
    hideydecorations!(relative_contributions_axis, ticks=false, ticklabels=false)
    hidexdecorations!(relative_contributions_axis)
    plot_raster!(raster_axis, raster)
    # rowsize!(sublayout, 1, Aspect)
    sublayout
end


fig3[1,3] = subfigC = let freq = 0.12, θ=-0, noise_amplitude = 0;
    raster_coords = IterTools.product(1:N_neurons, 1:N_times)
    raster_signal = map(raster_coords) do (i_neuron, i_time)
        (sin(2π*freq*i_time + θ*i_neuron/N_neurons) + 1) / 2 + (i_neuron % 10 == 0 ? 1 : 0)
    end
    raster_noise = noise_amplitude * rand(N_neurons, N_times)
    raster = (raster_signal + raster_noise) .> 0.5

   # contributions_axis = CairoMakie.Axis(fig3)
    relative_contributions_axis = CairoMakie.Axis(fig3)
    plot_relative_network_class_contributions!(relative_contributions_axis, raster, boundary, potential_contributions, lag_extents)
    raster_axis = CairoMakie.Axis(fig3)
    sublayout = GridLayout()
    title = Label(fig3, "phase shifting", tellheight=true, tellwidth=false)
    sublayout[:v] = [title, raster_axis, relative_contributions_axis]
    # plot_network_class_contributions!(contributions_axis, raster, lag_extents)
    # hideydecorations!(contributions_axis, ticks=false, ticklabels=false)
    # hidexdecorations!(contributions_axis)
    hideydecorations!(relative_contributions_axis, ticks=false, ticklabels=false)
    hidexdecorations!(relative_contributions_axis)
    plot_raster!(raster_axis, raster)
    # rowsize!(sublayout, 1, Aspect)
    sublayout
end

subfig3D = let freq = 0.12, θ=0, noise_amplitude=0.5, signal_amplitude=0.5;
    raster_coords = IterTools.product(1:N_neurons, 1:N_times)
    raster_signal = map(raster_coords) do (i_neuron, i_time)
        (sin(2π*freq*i_time + θ) + 1) / 2
    end
    raster_noise = rand(N_neurons, N_times)
    raster = (signal_amplitude .* raster_signal + noise_amplitude .* raster_noise) .> 0.5

   # contributions_axis = CairoMakie.Axis(fig3)
    relative_contributions_axis = CairoMakie.Axis(fig3)
    plot_relative_network_class_contributions!(relative_contributions_axis, raster, boundary, potential_contributions, lag_extents)
    raster_axis = CairoMakie.Axis(fig3)
    sublayout = GridLayout()
    title = Label(fig3, "SNR = 0dB", tellheight=true, tellwidth=false)
    sublayout[:v] = [title, raster_axis, relative_contributions_axis]

    # plot_network_class_contributions!(contributions_axis, raster, lag_extents)
    # hideydecorations!(contributions_axis, ticks=false, ticklabels=false)
    # hidexdecorations!(contributions_axis)
    plot_raster!(raster_axis, raster)
    sublayout
end

fig3[2,1] = subfig3D

subfig3E = let freq = 0.12, θ=0, noise_amplitude=2//3, signal_amplitude=1//3;
    raster_coords = IterTools.product(1:N_neurons, 1:N_times)
    raster_signal = map(raster_coords) do (i_neuron, i_time)
        (sin(2π*freq*i_time + θ) + 1) / 2
    end
    raster_noise = rand(N_neurons, N_times)
    raster = (signal_amplitude .* raster_signal + noise_amplitude .* raster_noise) .> 0.5

   # contributions_axis = CairoMakie.Axis(fig3)
    relative_contributions_axis = CairoMakie.Axis(fig3)
    plot_relative_network_class_contributions!(relative_contributions_axis, raster, boundary, potential_contributions, lag_extents)
    raster_axis = CairoMakie.Axis(fig3)
    sublayout = GridLayout()
    title = Label(fig3, "SNR = -6dB", tellheight=true, tellwidth=false)
    sublayout[:v] = [title, raster_axis, relative_contributions_axis]

    # plot_network_class_contributions!(contributions_axis, raster, lag_extents)
    # hideydecorations!(contributions_axis, ticks=false, ticklabels=false)
    # hidexdecorations!(contributions_axis)
    hideydecorations!(relative_contributions_axis, ticks=false, ticklabels=false)
    hidexdecorations!(relative_contributions_axis)
    plot_raster!(raster_axis, raster)
    sublayout
end

fig3[2,2] = subfig3E

subfig3F = let freq = 0.04, θ=0, noise_amplitude=1;
    raster_coords = IterTools.product(1:N_neurons, 1:N_times)
    raster_signal = map(raster_coords) do (i_neuron, i_time)
        (sin(2π*freq*i_time + θ) + 1) / 2
    end
    raster_noise = noise_amplitude * rand(N_neurons, N_times)
    raster = (raster_noise) .> 0.5

   # contributions_axis = CairoMakie.Axis(fig3)
    relative_contributions_axis = CairoMakie.Axis(fig3)
    plot_relative_network_class_contributions!(relative_contributions_axis, raster, boundary, potential_contributions, lag_extents)
    raster_axis = CairoMakie.Axis(fig3)
    sublayout = GridLayout()
    title = Label(fig3, "Noise", tellheight=true, tellwidth=false)
    sublayout[:v] = [title, raster_axis, relative_contributions_axis]

    # plot_network_class_contributions!(contributions_axis, raster, lag_extents)
    # hideydecorations!(contributions_axis, ticks=false, ticklabels=false)
    # hidexdecorations!(contributions_axis)
    hideydecorations!(relative_contributions_axis, ticks=false, ticklabels=false)
    hidexdecorations!(relative_contributions_axis)
    plot_raster!(raster_axis, raster)
    sublayout
end

fig3[2,3] = subfig3F

label_A = fig3[1,1,TopLeft()] = Label(fig3, "A", font=noto_sans_bold, textsize=56, halign=:left)
label_B = fig3[1,2,TopLeft()] = Label(fig3, "B", font=noto_sans_bold, textsize=56, halign=:left)
label_C = fig3[1,3,TopLeft()] = Label(fig3, "C", font=noto_sans_bold, textsize=56, halign=:left)
label_D = fig3[2,1,TopLeft()] = Label(fig3, "D", font=noto_sans_bold, textsize=56, halign=:left)
label_E = fig3[2,2,TopLeft()] = Label(fig3, "E", font=noto_sans_bold, textsize=56, halign=:left)
label_F = fig3[2,3,TopLeft()] = Label(fig3, "F", font=noto_sans_bold, textsize=56, halign=:left)


fig3

save(plotsdir("Fig3_$(boundary |> typeof)_$(neuron_max_lag)x$(time_max_lag)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")).png"), fig3)

end
end
