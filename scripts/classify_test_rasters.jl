include(srcdir("main.jl"))

using BenchmarkTools
using Test

@testset "Network Motif I" begin
    fn = datadir("I.csv")
    raster_I = load_raster(fn)

    contributions_I = triple_correlation_network_classifications(raster_I, 15, 15)
    @show contributions_I
    @test contributions_I[1] == 289 && all(contributions_I[begin+1:end] .== 0)
end

@testset "Network Motif II" begin
    fn = datadir("II.csv")
    raster_II = load_raster(fn)

    contributions_II = triple_correlation_network_classifications(raster_II, 15, 15)
    @show contributions_II
    @test contributions_II[2] > 0 && all(contributions_II[begin+2:end] .== 0)
end

@testset "Network Motif III" begin
    fn = datadir("III.csv")
    raster_III = load_raster(fn)

    contributions_III = triple_correlation_network_classifications(raster_III, 15, 15)
    @show contributions_III
    @test contributions_III[3] > 0 && all(contributions_III[begin+3:end] .== 0)
end

@testset "Network Motif IV" begin
    fn = datadir("IV.csv")
    raster_IV = load_raster(fn)

    contributions_IV = triple_correlation_network_classifications(raster_IV, 15, 15)
    @show contributions_IV
    @test contributions_IV[4] > 0 && all(contributions_IV[begin+4:end] .== 0)
end

@testset "Network Motif V" begin
    fn = datadir("V.csv")
    raster_V = load_raster(fn)

    contributions_V = triple_correlation_network_classifications(raster_V, 15, 15)
    @show contributions_V
    @test contributions_V[5] > 0 && all(contributions_V[begin+5:end] .== 0)
end