include(srcdir("main.jl"))

let fn = datadir("random_raster.csv")
    raster = load_raster(fn)

    triple_correlation_network_classifications(raster, 15, 15)
end