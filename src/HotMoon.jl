"""
Hot Moon Climate Simulator

A 1D/2D climate simulation for a moon using energy balance ODEs.

Uses Julia's multiple dispatch: create a moon with HotMoonBody(n_lat) for 1D
or HotMoonBody(n_lat, n_lon) for 2D, then call run_simulation(moon, hours, T0).
"""
module HotMoon

using ADTypes
using CoherentNoise
using Colors
using DifferentialEquations
using DifferentialEquations: PeriodicCallback
using Plots
using Printf
using Statistics

# =============================================================================
# Core Infrastructure
# =============================================================================

include("constants.jl")
include("performance.jl")
include("types.jl")
include("geometry.jl")
include("state.jl")
include("grid_iteration.jl")

# =============================================================================
# Physics Modules
# =============================================================================

include("physics/physics.jl")

# =============================================================================
# Solvers
# =============================================================================

include("solvers/solver_1d.jl")
include("solvers/solver_2d.jl")

# =============================================================================
# Visualization
# =============================================================================

include("visualization/visualization.jl")

# =============================================================================
# Output Management
# =============================================================================

include("output.jl")

# =============================================================================
# Exports: Types
# =============================================================================

export AbstractMoonBody, MoonBody1D, MoonBody2D
export HotMoonBody  # Convenience constructor

# =============================================================================
# Exports: Simulation Functions
# =============================================================================

export run_simulation
export run_simulation_with_moisture
export run_simulation_with_twolayer_atmosphere
export make_progress_callback


# =============================================================================
# Exports: Physics Functions (new names)
# =============================================================================

export compute_surface_albedo, compute_ice_fraction_from_temperature
export compute_total_ir_optical_depth, compute_ir_transmissivity_from_optical_depth
export compute_heat_transport_coefficient, compute_convection_activation
export compute_biome_blended_heat_capacity, classify_cell_biome
export compute_latitude_based_heat_capacity, compute_elevation_based_heat_capacity
export compute_saturation_moisture_at_temperature
export compute_ocean_evaporation_rate, compute_orographic_precipitation_rate
export compute_convective_ascent_rate, compute_total_descent_rate
export compute_ascent_moisture_transfer, compute_descent_moisture_transfer
export compute_descent_saturation_multiplier, compute_upper_mass_floor_restoration
export estimate_initial_temperature_and_moisture

# =============================================================================
# Exports: Terrain Functions
# =============================================================================

export create_fractal_noise_sampler, sample_terrain_elevation
export populate_heat_transport_coefficients!, populate_moisture_transport_coefficients!
export compute_directional_heat_transport_modifier, compute_directional_moisture_transport_modifier

# =============================================================================
# Exports: Biome Constants
# =============================================================================

export BIOME_OCEAN, BIOME_ICE_SHEET, BIOME_TUNDRA, BIOME_BOREAL_FOREST
export BIOME_COLD_STEPPE, BIOME_TEMPERATE_FOREST, BIOME_GRASSLAND
export BIOME_TROPICAL_FOREST, BIOME_SAVANNA, BIOME_HOT_DESERT
export BIOME_MOUNTAIN, BIOME_WETLAND, NUM_BIOMES, BIOME_NAMES, BIOME_HEAT_CAPACITY

# =============================================================================
# Exports: Geometry
# =============================================================================

export get_zenith, get_solar, is_eclipsed
export get_zenith_2d, get_solar_2d

# =============================================================================
# Exports: Visualization
# =============================================================================

# Variable types for dispatch
export PlotVariable, Temperature, Moisture, Precipitation
export UpperMass, UpperMoisture, Biome

# Generic 2D visualization
export plot_global_mean_timeseries, plot_field_snapshot
export plot_longitude_time_hovmoeller, plot_latitude_time_hovmoeller
export plot_latitude_mean_range, plot_elevation_map

# Biome visualization
export compute_biome_map, print_biome_statistics, plot_biome_with_legend

# Animation
export animate_field_evolution, animate_combined_climate, animate_all_variables

# 1D-specific visualization
export plot_1d_global_mean_temperature, plot_1d_latitude_temperature_profile
export plot_1d_temperature_heatmap, plot_1d_latitude_temperature_range
export plot_1d_latitude_temperature_timeseries, plot_1d_summary_to_files

# =============================================================================
# Exports: Performance
# =============================================================================

export set_threading, get_threading_status

# =============================================================================
# Exports: State Management
# =============================================================================

export unpack_temperature_field, unpack_moisture_field
export unpack_upper_mass_field, unpack_upper_moisture_field
export pack_temperature_moisture_state, pack_twolayer_state
export is_temperature_only_state, is_temperature_moisture_state, is_twolayer_state

# =============================================================================
# Exports: Grid Iteration
# =============================================================================

export foreach_cell, foreach_cell_sequential, foreach_cell_indexed
export foreach_valid_neighbor
export DIRECTION_NORTH, DIRECTION_SOUTH, DIRECTION_EAST, DIRECTION_WEST

# =============================================================================
# Exports: Output Management
# =============================================================================

export RunConfig, RunContext, RunInfo
export initialize_run, save_results!, generate_plots!, generate_animations!, finalize_run!
export load_run, list_runs, print_runs

# =============================================================================
# Exports: Constants
# =============================================================================

export ROTATION_PERIOD, ORBITAL_PERIOD, ECLIPSE_DURATION
export STEFAN_BOLTZMANN, SOLAR_CONSTANT
export BASE_IR_OPTICAL_DEPTH, WV_OPTICAL_SCALE, WV_REF_MOISTURE
export MOISTURE_REF_TEMP, MOISTURE_REF_SATURATION, CLAUSIUS_CLAPEYRON_SCALE
export EVAP_RATE, EVAP_THRESHOLD, PRECIP_RATE
export MOISTURE_DIFFUSION, MOISTURE_BARRIER_STRENGTH, ELEVATION_SCALE, LAPSE_RATE
export ELEVATION_GREENHOUSE_REDUCTION
export THERMAL_RESPONSE_SCALE, MIN_MOISTURE_CONVECTION, ASCENT_RATE_MAX
export BASE_DESCENT_RATE, MASS_DESCENT_COEFF, POLAR_SINK_STRENGTH, POLAR_SINK_LAT
export LIFT_TEMPERATURE_DROP, MOISTURE_SURVIVE_FRACTION, LIFT_FRACTION
export UPPER_MERIDIONAL_COEFF, UPPER_ZONAL_COEFF, WESTERLY_BIAS_STRENGTH
export DESCENT_DRYING_SCALE, U_FLOOR, U_INITIAL, M_UP_INITIAL
export DEFAULT_MOISTURE_ESTIMATE_KG_M2, INIT_MOISTURE_SATURATION_CAP

end # module
