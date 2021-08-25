
abstract type AbstractLag{N,T} end
function Base.Tuple(lag::AbstractLag)
    lag.λ
end
struct Lag{N,T} <: AbstractLag{N,T}
    λ::NTuple{N,T}
end
struct TernarizedLag{N,T} <: AbstractLag{N,T}# come up with better name
    λ::NTuple{N,T}
end
function ternarize_triplet_1d_lags((l1, l2))
    if l1 > l2 > 0
        return (2, 1)
    elseif l2 > l1 > 0
        return (1, 2)
    elseif l1 < l2 < 0
        return (-2, -1)
    elseif l2 < l1 < 0
        return (-1, -2)
    else
        return sign.((l1, l2))
    end
end

abstract type AbstractLagPair{N,T} end
function Base.Tuple(lag_pair::AbstractLagPair)
    tuple(Tuple(lag_pair.λ₁)..., Tuple(lag_pair.λ₂)...)
end
struct LagTriplet{N_dims,T_lag} <: AbstractLagPair{N_dims, T_lag}
    λ₁::Lag{N_dims,T_lag}
    λ₂::Lag{N_dims,T_lag}
end
struct LagMotif{N_dims,T_lag}  <: AbstractLagPair{N_dims, T_lag}
    λ₁::TernarizedLag{N_dims,T_lag}
    λ₂::TernarizedLag{N_dims,T_lag}
end
function LagMotif(lag_triplet::LagTriplet{N_dims,T_lag}) where {N_dims,T_lag}
    lag_pairs = zip(lag_triplet.λ₁.λ, lag_triplet.λ₂.λ)
    ternarized_1d_lag_pairs = ternarize_triplet_1d_lags.(lag_pairs)
    return LagMotif{N_dims,T_lag}(TernarizedLag.(zip(ternarized_1d_lag_pairs...))...)
end

struct LagMotif_TripleCorrelations{N_dims,T_val,T_lag,DCT<:AbstractDict{LagMotif{N_dims,T_lag},T_val}}
    prevalences::DCT
end
function reduce_to_lag_motifs(tc::TripleCorrelation{2,T_lag,4,T_val}) where {T_lag, T_val}
    prevalences = DefaultDict{LagMotif{2,T_lag},T_val}(zero(T_val))
    for coord ∈ IterTools.product(axes(tc)...)
        λ1 = Lag(coord[1:2])
        λ2 = Lag(coord[3:4])
        motif = LagMotif(LagTriplet(λ1, λ2))
        prevalences[motif] += tc[coord...]
    end
    return LagMotif_TripleCorrelations(prevalences)
end
Base.getindex(lm::LagMotif_TripleCorrelations, key) = getindex(lm.prevalences, key)
Base.setindex!(lm::LagMotif_TripleCorrelations, val, key) = setindex!(lm.prevalences, val, key)

