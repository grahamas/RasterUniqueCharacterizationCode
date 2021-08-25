

include(srcdir("main.jl"))

using BenchmarkTools


let raster_size = (150, 150),
        λ_max = (7, 7);

raster = generate_random_raster(raster_size)
plain_array = Array{Float64,4}(undef, length.([-λ_max[1]:λ_max[1], -λ_max[2]:λ_max[2], -λ_max[1]:λ_max[1], -λ_max[2]:λ_max[2]])...)
offset_array = OffsetArray{Float64}(copy(plain_array), -λ_max[1]:λ_max[1], -λ_max[2]:λ_max[2], -λ_max[1]:λ_max[1], -λ_max[2]:λ_max[2])

@btime _calculate_unscaled_triple_correlation!($plain_array, $raster, $λ_max)
@btime _calculate_unscaled_triple_correlation!($offset_array, $raster, $λ_max)

end