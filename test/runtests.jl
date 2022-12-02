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
end
