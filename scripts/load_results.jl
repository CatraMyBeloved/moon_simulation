"""
Load and analyze saved simulation results using the unified output system.

Usage:
    julia scripts/load_results.jl                  # Load latest run
    julia scripts/load_results.jl output/latest    # Load specific path
    julia scripts/load_results.jl --list           # List all runs
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Statistics

include(joinpath(@__DIR__, "../src/HotMoon.jl"))
using .HotMoon

# Handle command line arguments
if length(ARGS) > 0 && ARGS[1] == "--list"
    print_runs()
    exit(0)
end

# Determine which run to load
run_path = if length(ARGS) > 0
    ARGS[1]
else
    "output/latest"
end

if !isdir(run_path)
    println("ERROR: Run directory not found: $run_path")
    println("\nAvailable runs:")
    print_runs()
    exit(1)
end

# Load the run using unified API
println("="^60)
println("LOADING SIMULATION")
println("="^60)
println("  Path: $run_path")

sol, moon, metadata = load_run(run_path)

if sol === nothing
    error("No solution.jls found in $run_path")
end
if moon === nothing
    error("No moon.jls found in $run_path")
end

println("\nData loaded successfully!")
println("  Simulation time: $(sol.t[end]/3600) hours")
println("  Timesteps: $(length(sol.t))")

# Display statistics based on model type
println("\n" * "="^60)
println("ANALYSIS")
println("="^60)

if moon isa MoonBody1D
    println("  Model: 1D ($(moon.n_lat) latitude bands)")

    final_temps = sol.u[end] .- 273.15
    global_mean = sum(final_temps .* moon.cell_areas)

    println("\nFinal temperature statistics:")
    println("  Global mean: $(round(global_mean, digits=1))°C")
    println("  Equator: $(round(final_temps[1], digits=1))°C")
    println("  Pole: $(round(final_temps[end], digits=1))°C")

else
    # MoonBody2D - determine state type from solution size
    n_cells = moon.n_lat * moon.n_lon
    state_size = length(sol.u[end])

    if state_size == n_cells
        println("  Model: 2D temperature-only ($(moon.n_lat) × $(moon.n_lon))")
        T_final = reshape(sol.u[end], moon.n_lat, moon.n_lon)
        M_final = nothing

    elseif state_size == 2 * n_cells
        println("  Model: 2D with moisture ($(moon.n_lat) × $(moon.n_lon))")
        T_final = reshape(sol.u[end][1:n_cells], moon.n_lat, moon.n_lon)
        M_final = reshape(sol.u[end][n_cells+1:end], moon.n_lat, moon.n_lon)

    elseif state_size == 4 * n_cells
        println("  Model: 2D two-layer atmosphere ($(moon.n_lat) × $(moon.n_lon))")
        T_final = reshape(sol.u[end][1:n_cells], moon.n_lat, moon.n_lon)
        M_final = reshape(sol.u[end][n_cells+1:2n_cells], moon.n_lat, moon.n_lon)
        U_final = reshape(sol.u[end][2n_cells+1:3n_cells], moon.n_lat, moon.n_lon)
        M_up_final = reshape(sol.u[end][3n_cells+1:4n_cells], moon.n_lat, moon.n_lon)

        println("\nUpper layer statistics:")
        println("  Mean U: $(round(sum(U_final .* moon.cell_areas), digits=3))")
        println("  Mean M_up: $(round(sum(M_up_final .* moon.cell_areas), digits=4)) kg/m²")
    else
        error("Unknown state size: $state_size for grid $(moon.n_lat) × $(moon.n_lon)")
    end

    # Temperature stats
    T_final_C = T_final .- 273.15
    global_mean_T = sum(T_final_C .* moon.cell_areas)
    equator_idx = argmin(abs.(moon.latitudes))

    println("\nFinal temperature statistics:")
    println("  Global mean: $(round(global_mean_T, digits=1))°C")
    println("  Min temperature: $(round(minimum(T_final_C), digits=1))°C")
    println("  Max temperature: $(round(maximum(T_final_C), digits=1))°C")
    println("  Equator mean: $(round(mean(T_final_C[equator_idx, :]), digits=1))°C")

    # Moisture stats if available
    if M_final !== nothing
        global_mean_M = sum(M_final .* moon.cell_areas)
        println("\nFinal moisture statistics:")
        println("  Global mean: $(round(global_mean_M, digits=2)) kg/m²")
        println("  Min moisture: $(round(minimum(M_final), digits=2)) kg/m²")
        println("  Max moisture: $(round(maximum(M_final), digits=2)) kg/m²")
    end

    # Terrain info
    ocean_pct = 100 * sum(moon.elevation .< 0) / n_cells
    println("\nTerrain:")
    println("  Ocean coverage: $(round(ocean_pct, digits=1))%")
end

println("\n" * "="^60)
println("INTERACTIVE ANALYSIS")
println("="^60)
println("You can now interactively analyze the data:")
println("  - sol.t contains time points (seconds)")
println("  - sol.u contains state arrays at each time")
println("  - moon contains the spatial grid and terrain")
println("\nExample: unpack_temperature_field(sol.u[end], moon)")
println("="^60)
