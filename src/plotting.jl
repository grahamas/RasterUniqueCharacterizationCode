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

function plot_network_class_contributions!(ax, raster::BitMatrix, boundary, max_lags; fillto=10^-1)

    network_class_contributions = sequence_class_tricorr_unrolled(raster, boundary, max_lags)
    fillto=if ax.yscale[] == log10
        ax.ytickformat[] = Makie.automatic
        fillto_exp = minimum(floor.(log10.(network_class_contributions)))
        10^fillto_exp
    else
        0.
    end
    plt = barplot!(ax, 1:14 |> collect, network_class_contributions; fillto=fillto);
    ax.xticks = [1:5:14...]
    ax.xtickformat[] = xs -> (roman_encode ∘ Int).(xs)
    tightlimits!(ax)
    ax.ylabel[] = "number of motifs"
    ax.xlabel[] = "network motif"
    plt
end

using Makie
function plot_relative_network_class_contributions!(ax, raster::BitMatrix, boundary, potential_contributions, max_lags)
    ax.ytickformat[] = Makie.automatic
    network_class_contributions = sequence_class_tricorr_unrolled(raster, boundary, max_lags)
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
    #ax.tellheight = true
    ax.bottomspinevisible = false
    ax.leftspinevisible = false
    ax.topspinevisible = false
    ax.rightspinevisible = false
    ax.autolimitaspect = 1.
end

using Distributions
function rand_spike_train(λ::Number, n_bins)
    train = zeros(Int, n_bins)
    train .= 0
    rand_spike_train!(train, λ)
    return train
end
function rand_spike_train!(train, λ)
    n_bins = length(train)
    next = ceil(Int, rand(Exponential(1/λ)))
    while next <= n_bins
        train[next] = 1
        next += ceil(Int, rand(Exponential(1/λ)))
    end
end

function rand_raster(λ, n_channels, n_bins)
    raster = zeros(Int, n_channels, n_bins)
    for i_ch ∈ 1:n_channels
        @views rand_spike_train!(raster[i_ch, :], λ)
    end
    return raster
end
