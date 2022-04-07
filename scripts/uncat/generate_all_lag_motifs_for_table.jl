@quickactivate "RasterUniqueCharacterizationCode"

using DataFrames, CSV

include(srcdir("motif_classes.jl"))
include(srcdir("enumerate_lag_motifs_OLD.jl"))
include(srcdir("enumerate_lag_motifs.jl"))
include(srcdir("latex_table.jl"))

classes = generate_all_lag_motif_classes()
# remap_old_filenames(plotsdir("old_motif_figs"), plotsdir("motif_figs"))
info_flow_classes = info_flow_classify_lag_motif_classes(classes)
df = DataFrame(prettify_row.(name_info_flow_row.(info_flow_classes)))
CSV.write(datadir("info_flow_classified_table_abs.csv"), df)

open(datadir("latex_table.tex"), "w") do io
    write(io, latexify(df))
end