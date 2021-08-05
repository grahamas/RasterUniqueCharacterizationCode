using DataFrames, CSV


include(srcdir("motif_classes.jl"))
include(srcdir("latex_table.jl"))

classes = generate_all_motif_classes()
classified = connectivity_classify_motif_classes(classes)
df = DataFrame(prettify_row.(name_classified_row.(classified)))
CSV.write(datadir("classified_table_abs.csv"), df)

open(datadir("latex_table.tex"), "w") do io
    write(io, latexify(df))
end