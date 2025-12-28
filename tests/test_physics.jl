"""
Unit tests for physics functions
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test

include(joinpath(@__DIR__, "../src/HotMoon.jl"))
using .HotMoon

@testset "HotMoon Physics Tests" begin
    moon = HotMoonBody(18)

    @testset "Albedo" begin
        # Cold should have higher albedo than warm (ice-albedo feedback)
        albedo_cold = get_albedo(250.0)
        albedo_warm = get_albedo(300.0)
        @test albedo_cold >= albedo_warm

        # Albedo should be between 0 and 1
        @test 0.0 <= albedo_cold <= 1.0
        @test 0.0 <= albedo_warm <= 1.0
    end

    @testset "Greenhouse" begin
        # Warm should have stronger greenhouse than cold (water vapor)
        gh_cold = get_greenhouse(250.0)
        gh_warm = get_greenhouse(300.0)
        @test gh_warm >= gh_cold

        # Greenhouse should be between 0 and 1
        @test 0.0 <= gh_cold <= 1.0
        @test 0.0 <= gh_warm <= 1.0
    end

    @testset "Solar Geometry" begin
        # With 28-hour rotation: phase=0 is noon, phase=0.25 is sunset, phase=0.5 is midnight
        # Daytime: phase < 0.25 or phase > 0.75 → t < 7h or t > 21h

        # Daytime should have zenith angle (t=2h, phase=0.071, morning)
        t_day = 2.0 * 3600
        zenith = get_zenith(t_day, 0.0)
        @test zenith !== nothing

        # Nighttime should return nothing (t=14h, phase=0.5, midnight)
        t_night = 14.0 * 3600
        zenith_night = get_zenith(t_night, 0.0)
        @test zenith_night === nothing
    end

    @testset "Solar Heating" begin
        # Daytime should give positive heating (t=2h is morning)
        Q_day = get_solar(2.0*3600, 0.0, 285.0)
        @test Q_day > 0

        # Nighttime should give zero heating (t=14h is midnight)
        Q_night = get_solar(14.0*3600, 0.0, 285.0)
        @test Q_night == 0.0

        # Equator should get more than high latitude at same local time
        Q_equator = get_solar(2.0*3600, 0.0, 285.0)
        Q_highalt = get_solar(2.0*3600, 80.0, 285.0)
        @test Q_equator >= Q_highalt
    end

    @testset "Heat Capacity Zones" begin
        # Equatorial zone should have highest heat capacity
        C_equator = get_heat_capacity(15.0)  # 15° is in equatorial zone
        C_midlat = get_heat_capacity(45.0)   # 45° is in mid-latitude zone
        C_polar = get_heat_capacity(75.0)    # 75° is in polar zone

        @test C_equator > C_midlat
        @test C_midlat > C_polar

        # All should be positive
        @test C_equator > 0
        @test C_midlat > 0
        @test C_polar > 0

        # Test zone boundaries
        @test get_heat_capacity(29.0) == get_heat_capacity(15.0)  # Both equatorial
        @test get_heat_capacity(31.0) == get_heat_capacity(45.0)  # Both mid-lat
        @test get_heat_capacity(61.0) == get_heat_capacity(85.0)  # Both polar
    end
end

println("\nAll tests passed! ✓")