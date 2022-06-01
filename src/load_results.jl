@quickactivate "RasterUniqueCharacterizationCode"

using TripleCorrelations
using DataFrames
using JLD2

function load_results_df(filename)
    loaded = JLD2.load(filename)
    prior_results_dict = loaded["prior_results_dict"]
    DataFrame(
        mapreduce(vcat, keys(prior_results_dict)) do key
            _, _, peristimulus_results = prior_results_dict[key]
            mapreduce(vcat, peristimulus_results) do result
                map(axes(result.mean_effect,1)) do i_motif
                    (
                        signal_motif=key.motif_class,
                        detect_motif=TripleCorrelations.offset_motif_numeral(i_motif),
                        mean_effect=result.mean_effect[i_motif], 
                        #std_effect=result.std_effect[i_motif],
                        proportion_rejected=result.proportion_rejected[i_motif],
                        n_trials=result.sample_size,
                        n_hypothesis_tests=key.n_resamples
                    )
                end
            end
        end
    )
end
