include(srcdir("rasters.jl"))

let fn = datadir("random_raster.csv")
    raster = random_raster(150, 150)
    write_raster(fn, raster)
end