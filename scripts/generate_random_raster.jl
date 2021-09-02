include(srcdir("rasters.jl"))

let fn = datadir("random_raster_150x150.csv")
    raster = generate_random_raster((150, 150))
    write_raster(fn, raster)
end