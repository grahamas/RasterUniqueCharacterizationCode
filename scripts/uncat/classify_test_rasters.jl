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

@testset "Network Motif VI" begin
    fn = datadir("VI.csv")
    raster_VI = load_raster(fn)

    contributions_VI = triple_correlation_network_classifications(raster_VI, 15, 15)
    @show contributions_VI
    @test contributions_VI[6] > 0 && all(contributions_VI[begin+6:end] .== 0)
end

@testset "Network Motif VII" begin
    fn = datadir("VII.csv")
    raster_VII = load_raster(fn)

    contributions_VII = triple_correlation_network_classifications(raster_VII, 15, 15)
    @show contributions_VII
    @test contributions_VII[7] > 0 && all(contributions_VII[begin+7:end] .== 0)
end

@testset "Network Motif VIII" begin
    fn = datadir("VIII.csv")
    raster_VIII = load_raster(fn)

    contributions_VIII = triple_correlation_network_classifications(raster_VIII, 15, 15)
    @show contributions_VIII
    @test contributions_VIII[8] > 0 && all(contributions_VIII[begin+8:end] .== 0)
end

@testset "Network Motif IX" begin
    fn = datadir("IX.csv")
    raster_IX = load_raster(fn)

    contributions_IX = triple_correlation_network_classifications(raster_IX, 15, 15)
    @show contributions_IX
    @test contributions_IX[9] > 0 && all(contributions_IX[begin+9:end] .== 0)
end

@testset "Network Motif X" begin
    fn = datadir("X.csv")
    raster_X = load_raster(fn)

    contributions_X = triple_correlation_network_classifications(raster_X, 15, 15)
    @show contributions_X
    @test contributions_X[10] > 0 && all(contributions_X[begin+10:end] .== 0)
end

@testset "Network Motif XI" begin
    fn = datadir("XI.csv")
    raster_XI = load_raster(fn)

    contributions_XI = triple_correlation_network_classifications(raster_XI, 15, 15)
    @show contributions_XI
    @test contributions_XI[11] > 0 && all(contributions_XI[begin+11:end] .== 0)
end

@testset "Network Motif XII" begin
    fn = datadir("XII.csv")
    raster_XII = load_raster(fn)

    contributions_XII = triple_correlation_network_classifications(raster_XII, 15, 15)
    @show contributions_XII
    @test contributions_XII[12] > 0 && all(contributions_XII[begin+12:end] .== 0)
end

@testset "Network Motif XIII" begin
    fn = datadir("XIII.csv")
    raster_XIII = load_raster(fn)

    contributions_XIII = triple_correlation_network_classifications(raster_XIII, 15, 15)
    @show contributions_XIII
    @test contributions_XIII[13] > 0 && all(contributions_XIII[begin+13:end] .== 0)
end

@testset "Network Motif XIV" begin
    fn = datadir("XIV.csv")
    raster_XIV = load_raster(fn)

    contributions_XIV = triple_correlation_network_classifications(raster_XIV, 15, 15)
    @show contributions_XIV
    @test contributions_XIV[14] > 0 && all(contributions_XIV[begin+14:end] .== 0)
end