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

@load "/home/graham/git/RasterUniqueCharacterizationCode/remote_plots/plots/rate_divide_200trials_3sigs_IND_periodic_2022_05_17-124349/results.jld2"

motif_class = "I"
boundary_class = Periodic

prior_keys = keys(prior_results_dict) |> collect
mc_key = prior_keys[findfirst(x -> (x.motif_class == motif_class && x.boundary isa boundary_class), prior_keys)]

l_trials_epochs, trialavg_raster, peristimulus_results = prior_results_dict[mc_key]


peristimulus_results_by_motif = DataFrame(mapreduce(vcat, peristimulus_results) do result
    map(axes(result.mean_effect,1)) do i_motif
        (
            signal_motif=motif_class,
            detect_motif=TripleCorrelations.offset_motif_numeral(i_motif),
            mean_effect=result.mean_effect[i_motif], 
            #std_effect=result.std_effect[i_motif],
            proportion_rejected=result.proportion_rejected[i_motif],
            sample_size=result.sample_size
        )
    end
end)

plt_power = data(peristimulus_results_by_motif) * mapping(:sample_size, :proportion_rejected, color=:detect_motif) * (visual(Scatter) + smooth(span=0.9, degree=2))
color_list = distinguishable_colors(N_MOTIFS, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
f_power = draw(plt_power, palettes=(color=color_list,))#, axis=(; title="Signal motif-class $(motif_class)"))