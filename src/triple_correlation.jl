# FIXME should encode stepsize
struct TripleCorrelation{N_dims,T_lag,N_dims2,T_val,ARR<:AbstractArray} <: AbstractArray{T_val,N_dims2}
    arr::ARR
    extent::NTuple{N_dims,T_lag}
    function TripleCorrelation(arr::ARR, extent::NTuple{N_dims,T_lag}) where {N_dims,N_dims2,T_val,T_lag,ARR<:AbstractArray{T_val,N_dims2}}
        N_dims * 2 == N_dims2 || throw(ArgumentError("Triple correlation should have twice as many dimensions as source."))
        new{N_dims,T_lag,N_dims2,T_val,ARR}(arr, extent)
    end
end
Base.parent(tc::TripleCorrelation) = tc.arr
Base.size(tc::TripleCorrelation) = size(parent(tc))
Base.getindex(tc::TripleCorrelation, i::Int) = getindex(parent(tc), i)
Base.getindex(tc::TripleCorrelation, I::Vararg{Int, N}) where N = getindex(parent(tc), I...)
Base.setindex!(tc::TripleCorrelation, v, i::Int) = setindex!(parent(tc), v, i)
Base.setindex!(tc::TripleCorrelation, v, I::Vararg{Int, N}) where N = setindex!(parent(tc), v, I...)
Base.IndexStyle(::Type{TripleCorrelation{N,T,N2,Tv,ARR}}) where {N,T,N2,Tv,ARR}= Base.IndexStyle(ARR)
Base.axes(tc::TripleCorrelation) = axes(parent(tc))

function triple_correlation(raster::R, λ_max::NTuple{N,T_lag}) where {T_src,N,T_lag,R<:AbstractArray{T_src,N}}
    TripleCorrelation(_calculate_scaled_triple_correlation(raster, λ_max), size(raster))
end



function _calculate_unscaled_triple_correlation!(correlation::Array{T_cor,4}, raster::R, λ_max::NTuple{2,T_lag}) where {T_cor,T_src,T_lag,R<:AbstractArray{T_src,2}}
    single_λ_ranges = [-l:l for l ∈ λ_max]
    neuron_lag_range, time_lag_range = single_λ_ranges


    (N_neurons, N_times) = size(raster)
    time_range = (1-minimum(time_lag_range)):(N_times-maximum(time_lag_range))
    neuron_range = (1-minimum(neuron_lag_range)):(N_neurons-maximum(neuron_lag_range))

    @turbo for i_n1 ∈ eachindex(neuron_lag_range), i_t1 ∈ eachindex(neuron_lag_range),  i_n2 ∈ eachindex(time_lag_range), i_t2 ∈ eachindex(time_lag_range)
        n1, n2, t1, t2 = neuron_lag_range[i_n1], neuron_lag_range[i_n2], 
time_lag_range[i_t1], time_lag_range[i_t2]
        accum = 0
        for i_neuron ∈ neuron_range, i_time ∈ time_range
            accum += raster[i_neuron, i_time] * raster[i_neuron+n1,i_time+t1] * raster[i_neuron+n2,i_time+t2]
        end
        correlation[i_n1,i_t1,i_n2,i_t2] = accum
    end
end

function _calculate_unscaled_triple_correlation!(correlation::OffsetArray{T_cor,4}, raster::R, λ_max::NTuple{2,T_lag}) where {T_cor,T_src,T_lag,R<:AbstractArray{T_src,2}}
    _calculate_unscaled_triple_correlation!(parent(correlation), raster, λ_max)
end


function _calculate_unscaled_triple_correlation(raster::BitMatrix, λ_max::Tuple)
    λ_ranges = Tuple(repshoweat([-l:l for l ∈ λ_max], outer=(2,)))
    correlation = zeros(Float64, λ_ranges)

    _calculate_unscaled_triple_correlation!(correlation, raster, λ_max)
    return correlation
end

function _calculate_scaled_triple_correlation(raster::BitMatrix, λ_max::Tuple)
    unscaled_correlation = _calculate_unscaled_triple_correlation(raster, λ_max)
    unscaled_correlation ./= calculate_scaling_factor(raster, λ_max)
end
_calculate_scaled_triple_correlation(raster::OffsetArray, args...) = _calculate_scaled_triple_correlation(parent(raster), args...)

_hi_bound(l1, l2, M) = M - max(0, l1, l2)
_lo_bound(l1, l2) = 1 + min(0, l1, l2)

function n_triplets(n1, t1, n2, t2, N, T)
    (_hi_bound(n1, n2, N) - _lo_bound(n1, n2) + 1) *
    (_hi_bound(t1, t2, T) - _lo_bound(t1, t2) + 1)
end


####### High dimensional fallback #######


function _calculate_unscaled_triple_correlation!(correlation::OffsetArray{T_cor,N2}, raster::R, λ_max::NTuple{N,T_lag}) where {T_cor,N2,T_src,N,T_lag,R<:AbstractArray{T_src,N}}
    @warn "Using slow fallback for >2D image"
    N * 2 == N2 || throw(ArgumentError("Triple correlation should have twice as many dimensions as source."))
    λ_ranges = [-l:l for l ∈ λ_max]
    λs = product(λ_ranges...)
    inner_raster_coords = CartesianIndex(1 .+ λ_max):CartesianIndex(size(raster) .- λ_max)
    for i_λ₁ ∈ eachindex(λs), i_λ₂ ∈ eachindex(λs)
        λ₁ = λs[i_λ₁]
        λ₂ = λs[i_λ₂]
        accum = 0
        for i_base_coord ∈ eachindex(inner_raster_coords)
            base_coord = inner_raster_coords[i_base_coord]
            accum += raster[base_coord] * 
                     raster[CartesianIndex(Tuple(base_coord) .+ λ₁)] *
                     raster[CartesianIndex(Tuple(base_coord) .+ λ₂)]
        end
        correlation[λ₁..., λ₂...] = accum
    end
end
