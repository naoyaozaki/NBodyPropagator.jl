using NBodyProblem
import SPICE
using Test

@testset "NBodyProblem.jl" begin

    @testset "Case 1: Test values of AU,ID,GM,RE,NAME" begin
        # Compared with nomad (Python) results
        @test ssd.AU == 149597870.7 # km
        @test ssd.ID["EARTH"] == 399 # <Integer>
        @test ssd.ID["PHOBOS"] == 401 # <Integer>
        @test ssd.GM["EARTH"] == 3.98600435436095925979e+05 # km3/s2
        @test ssd.GM["MOON"] == 4.90280006616379614570e+03 # km3/s2
        @test ssd.GM["TITAN"] == 8.97813884530737595924e+03 # km3/s2
        @test ssd.RE["EARTH"] == 6.37813659999999981665e+03 # km
        @test ssd.RE["VENUS"] == 6.05180000000000018190e+03 # km
        @test ssd.RE["TITAN"] == 2.57515000000000009095e+03 # km
        @test ssd.NAME[399] == "EARTH" # <String>
        @test ssd.NAME[501] == "IO" # <String>
    end

    @testset "Case 2: Test values of OE (Orbital Elements)" begin
        # Compared with nomad (Python) results
        @test ssd.OE["JUPITER"]["sma"] == 7.78340816692710757256e+08 # km
        @test ssd.OE["JUPITER"]["ecc"] == 4.83862399999999970301e-02 # -
        @test ssd.OE["JUPITER"]["inc"] == (π / 180) * 1.30439695000000011049e+00 # rad
        @test ssd.OE["JUPITER"]["meanlong"] == (π / 180) * 3.43964405099999979143e+01 # rad
        @test ssd.OE["JUPITER"]["lnode"] == (π / 180) * 1.00473909090000006472e+02 # rad
        @test ssd.OE["JUPITER"]["argp"] == (π / 180) * 1.47284798299999994953e+01 # rad
    end

    @testset "Case 3: Test values of SPICE" begin
        # Compared with jTOP (MATLAB) results
        init_spice_kernels()

        # Moon state vector
        et = SPICE.str2et("2022/12/02 13:52:32 UTC") # Initial Epoch
        expected = [5.07170845744684636593e+07, 1.38670086263312786818e+08, -2.87103076627328991890e+04, -2.85197938092883269690e+01, 1.10662945701793731246e+01, 6.97827903405290328465e-02]
        @test SPICE.spkez(301, et, "ECLIPJ2000", "NONE", 10)[1] == expected

        # Europa state vector
        et = SPICE.str2et("2025/12/02 13:52:32 TDB") # Initial Epoch
        expected = [-6.71573226872288389131e+05, 4.57777683529129062663e+04, 5.90781640300525486964e+03, -1.01543406126799906630e+00, -1.23500573487305977949e+01, -5.84209720106468921585e+00]
        @test SPICE.spkez(502, et, "J2000", "NONE", 599)[1] == expected
    end
end
