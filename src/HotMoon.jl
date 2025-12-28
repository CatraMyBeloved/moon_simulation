"""
Hot Moon Climate Simulator

A 1D/2D climate simulation for a tidally-locked moon using energy balance ODEs.

Uses Julia's multiple dispatch: create a moon with HotMoonBody(n_lat) for 1D
or HotMoonBody(n_lat, n_lon) for 2D, then call run_simulation(moon, hours, T0).
"""
module HotMoon

using CoherentNoise
using Colors
using DifferentialEquations
using Plots
using Statistics

# Include source files in dependency order
include("constants.jl")
include("types.jl")
include("physics.jl")
include("geometry.jl")
include("solver_1d.jl")
include("solver_2d.jl")
include("visualization.jl")

# Export types
export AbstractMoonBody, MoonBody1D, MoonBody2D
export HotMoonBody  # Convenience constructor

# Export unified simulation function (dispatches on moon type)
export run_simulation

# Export physics functions
export get_albedo, get_greenhouse, get_heat_capacity, get_transport_coefficient
export generate_elevation, create_elevation_sampler
export calculate_directional_transport, calculate_transport_coefficients!

# Export geometry functions
export get_zenith, get_solar, is_eclipsed
export get_zenith_2d, get_solar_2d

# Export visualization - 1D
export plot_global_mean, plot_latitude_profile, plot_heatmap
export plot_latitude_mean_range, plot_latitude_timeseries, plot_summary

# Export visualization - 2D
export plot_longitude_timeseries, plot_latitude_mean_range_2d
export plot_temperature_heatmap_2d, plot_hovmoeller, plot_polar_projection
export plot_elevation_map

# Export constants that users might want to reference
export ROTATION_PERIOD, ORBITAL_PERIOD, ECLIPSE_DURATION
export STEFAN_BOLTZMANN, SOLAR_CONSTANT

end # module
