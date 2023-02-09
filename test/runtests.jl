using NBodyPropagator
import SPICE
using Test

@testset "NBodyPropagator.jl" begin

    # Get Solar System Dynamics Constant
    ssd = SolarSystemDynamics()

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
        # Moon state vector
        et = SPICE.str2et("2022/12/02 13:52:32 UTC") # Initial Epoch
        expected = [5.07170845744684636593e+07, 1.38670086263312786818e+08, -2.87103076627328991890e+04, -2.85197938092883269690e+01, 1.10662945701793731246e+01, 6.97827903405290328465e-02]
        @test SPICE.spkez(301, et, "ECLIPJ2000", "NONE", 10)[1] == expected

        # Europa state vector
        et = SPICE.str2et("2025/12/02 13:52:32 TDB") # Initial Epoch
        expected = [-6.71573226872288389131e+05, 4.57777683529129062663e+04, 5.90781640300525486964e+03, -1.01543406126799906630e+00, -1.23500573487305977949e+01, -5.84209720106468921585e+00]
        @test SPICE.spkez(502, et, "J2000", "NONE", 599)[1] == expected
    end

    @testset "Case 4: Interplanetary N-body propagation" begin
        # Expected solution computed by jTOP
        expected = [+4.223450194979922e+07, +1.167640417577651e+08, +7.759087742359675e+06, +9.233242342590227e+00, -5.023717927929555e+00, +2.053384118051583e+00]  # km, km/s

        # Parameter Setting
        list_bodies = [10, 399, 301, 299, 499, 599]

        # Initial Condition
        et0 = SPICE.str2et("2019 JAN 01 12:00:00 TDB")
        x0 = [1.0e7, 1.0e8, 1.0e6, 15.0, 20.0, 3.0]
        tspan = (et0, et0 + 30.0 * 86400.0)

        # Integrate
        nbp = NBodyProblem(x0, tspan, list_bodies)
        sol = propagate(nbp)

        @test sol[:, end] ≈ expected rtol = 1e-10
    end

    @testset "Case 5: Earth orbiting N-body propagation" begin
        # Expected solution computed by jTOP
        expected = [-2.291792316114288e+04, +4.609531317147258e+04, +8.960615353005227e+03, -2.217738317104420e+00, -1.185657482487868e+00, -2.299275473222509e-01]  # km, km/s

        # Parameter Setting
        list_bodies = [10, 399, 301, 299, 499, 599]

        # Initial Condition
        et0 = SPICE.str2et("2023/07/01 00:00:00 UTC")
        x0 = [42000, 10, 50, 0.5, 3.080663355435613, 0.6]
        tspan = (et0, et0 + 30.0 * 86400.0)

        # Integrate
        nbp = NBodyProblem(x0, tspan, list_bodies, id_center=399, ref_frame="J2000", lsf=1.0e5, tsf=1.0e5)
        sol = propagate(nbp)

        @test sol[:, end] ≈ expected rtol = 1e-6
    end

    @testset "Case 6: Interplanetary N-body propagation (SSB center) with STM" begin

        # Expected solution computed by jTOP
        expected_sv = [5.775330136689539e+07, 7.778671850845341e+07, -3.473805401858905e+03, -3.064706431615033e+01, 1.995206280152122e+01, -2.462928190848737e-03]  # km, km/s
        expected_stm = [
            1.847376983668235e+00 3.273089387090845e-01 1.051837686670713e-04 3.185856815602886e+06 4.319076300587049e+05 7.913848923152865e+01
            4.320189665583292e-01 7.408999038578260e-01 3.527313765385191e-05 4.796308155139408e+05 2.475618183196456e+06 3.830410284418312e+01
            1.041336022608508e-04 2.632574832527860e-05 5.780670629196961e-01 7.865927745045660e+01 3.422589891588510e+01 2.209053865236973e+06
            5.987311543450867e-07 3.341269682700137e-07 7.396314793709315e-11 1.545103188666622e+00 6.030866314541509e-01 7.809866498497458e-05
            5.321348585812730e-07 -6.256611273539320e-08 4.302221806248853e-11 7.136537012085610e-01 1.069006124741629e+00 5.587006133725106e-05
            7.197456461855974e-11 2.610223026701529e-11 -3.060610610748985e-07 7.698701383399397e-05 4.642107311960450e-05 5.603063175989961e-01
        ] # km, km/s
        expected_dxdt0 = [1.056147717588803e-02, 6.779051064188787e-03, 1.495088355184512e-04, 7.177542471178463e-09, 7.531287238037084e-09, 1.089708655901337e-10]

        # Parameter Setting
        list_bodies = [10, 399, 301, 299, 499, 599]

        # Initial Condition
        et0 = SPICE.str2et("2019 JAN 01 12:00:00 TDB")
        x0 = [1.0e8, 0.0, 0.0, 0.0, 35.0, 0.0]
        tspan = (et0, et0 + 30.0 * 86400.0)

        # Integrate
        nbp = NBodyProblem(x0, tspan, list_bodies, id_center=0, need_stm=true)
        state_all, stm_all, dxdt0_all = propagate(nbp)

        @test state_all[:, end] ≈ expected_sv rtol = 1e-5
        @test stm_all[:, :, end] ≈ expected_stm rtol = 1e-6
        @test dxdt0_all[:, end] ≈ expected_dxdt0 rtol = 1e-5

    end

    @testset "Case 7: Earth-orbiting N-body propagation (Earth center) with STM" begin
        # Expected solution computed by jTOP
        expected_sv = [-2.076802656535074e+03, 4.195026441409034e+04, -3.724592780590125e+00, -3.076821569056789e+00, -1.522475993300681e-01, -6.201698518781910e-04]  # k km, km/s
        expected_stm = [
            +5.676027474516669e+02 +1.019340915416001e+00 +2.714353988420178e-04 +2.774302596792756e+04 +7.711356579263460e+06 -2.235969046951308e+01
            +3.013664988602837e+01 +1.056255718751322e+00 +4.859255490766157e-06 +1.508648376854883e+04 +4.101833843734498e+05 +5.174580428182746e-02
            +1.138431767051646e-01 +1.250051853981020e-04 -5.245703996226803e-02 +4.480636984153578e+00 +1.546659652345703e+03 +1.359941126640226e+04
            -1.988671974500659e-03 +6.986859990077905e-05 +1.374723705414486e-08 +9.023001865252026e-01 -2.701564879310958e+01 +8.086374067414333e-05
            +4.171330552995736e-02 +7.881340971863755e-05 +1.331341475367478e-08 +2.090054618527219e+00 +5.677108329588499e+02 -1.439267396309650e-03
            -3.688142204328751e-06 +8.361081474993552e-09 -7.332738614462241e-05 +2.636725846778809e-05 -5.018640526647297e-02 -5.319394856551071e-02
        ]  # km, km/s
        expected_dxdt0 = [1.113650997738780e-03, 5.941571232254622e-05, 2.905738548207350e-06, -3.920574819395139e-09, 8.173154187726920e-08, -9.277491401613935e-11]

        # Parameter Setting
        list_bodies = [399, 301, 10, 499, 599]

        # Initial Condition
        et0 = SPICE.str2et("2016 NOV 01 00:00:00 UTC")
        x0 = [42000.0, 0.0, 0.0, 0.0, sqrt(398600.4418 / 42000.0), 0.0]
        tspan = (et0, et0 + 30.0 * 86400.0)

        # Integrate
        nbp = NBodyProblem(x0, tspan, list_bodies, id_center=399, need_stm=true, lsf=1.0, tsf=1.0, msf=1.0)
        state_all, stm_all, dxdt0_all = propagate(nbp)

        @test state_all[:, end] ≈ expected_sv rtol = 1e-5
        @test stm_all[:, :, end] ≈ expected_stm rtol = 1e-6
        @test dxdt0_all[:, end] ≈ expected_dxdt0 rtol = 1e-5

    end

end
