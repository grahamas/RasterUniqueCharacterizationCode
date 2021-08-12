noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")

abbrev_count_label = x -> begin
        if x >= 1000000000
            try
                "$(Int(x / 1000000000))B"
            catch
                "$(x / 1000000000)B"
            end
        elseif x >= 1000000
            try
                "$(Int(x / 1000000))M"
            catch
                "$(x / 1000000)M"
            end
        elseif x >= 1000
            try
                "$(Int(x / 1000))K"
            catch
                "$(x / 1000)K"
            end
        else
            "$(Int(x))"
        end
    end

bar_theme = Theme(
    fontsize=30,
    Axis = (
        backgroundcolor = :white,
        leftspinevisible = true,
        rightspinevisible = false,
        bottomspinevisible = true,
        topspinevisible = false,
        xgridcolor = :white,
        ygridcolor = :white,
        ytickformat = xs -> abbrev_count_label.(xs)
    )
)

function plot_network_class_contributions!(ax, raster::BitMatrix, neuron_window_size, time_window_size; fillto=0.)

    network_class_contributions = triple_correlation_network_classifications(raster, neuron_window_size, time_window_size)
    plt = barplot!(ax, 1:14 |> collect, network_class_contributions; fillto=fillto);
    ax.xticks = [1:5:14...]
    ax.xtickformat[] = xs -> (roman_encode ∘ Int).(xs)
    tightlimits!(ax)
    ax.ylabel[] = "number of motifs"
    ax.xlabel[] = "network motif"
    plt
end

using Makie
function plot_relative_network_class_contributions!(ax, raster::BitMatrix, potential_contributions, neuron_window_size, time_window_size)
    ax.ytickformat[] = Makie.automatic
    network_class_contributions = triple_correlation_network_classifications(raster, neuron_window_size, time_window_size)
    relative_contributions = network_class_contributions ./ potential_contributions
    plt = stem!(ax, 1:14 |> collect, (relative_contributions); ytickformat=identity);
    ax.xticks = [1:5:14...]
    ax.xtickformat[] = xs -> (roman_encode ∘ Int).(xs)
    tightlimits!(ax)
    ax.ylabel[] = "proportion of potential motifs"
    ax.xlabel[] = "network motif"
    xlims!(ax, (-0.5, 14.5))
    ylims!(ax, (0., 1.))
    plt
end

function plot_raster!(ax, raster::BitMatrix)
    plt = heatmap!(ax, raster', colormap=:binary)
    hidedecorations!(ax)
    ax.aspect = DataAspect()
    ax.tellheight = true
    ax.bottomspinevisible = false
    ax.leftspinevisible = false
    plt
end