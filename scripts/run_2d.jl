"""
Main script to run 2D Hot Moon simulation
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Statistics

include(joinpath(@__DIR__, "../src/HotMoon.jl"))
using .HotMoon

println("="^60)
println("Hot Moon 2D Climate Simulation")
println("="^60)

# Create moon with 90 lat × 180 lon grid (2° × 2° resolution)
println("\nCreating moon...")
moon = HotMoonBody(90, 180, seed=32, sea_level=0.0, scale=0.01, octaves=5)
println(moon)
println("  Total cells: $(moon.n_lat * moon.n_lon)")

# Initial conditions: uniform 285K (12°C)
T0 = fill(285.0, moon.n_lat, moon.n_lon)
println("Initial temperature: $(T0[1,1] - 273.15)°C (uniform)")

# Run simulation
sim_hours = 10000.0
println("\nRunning simulation for $(round(Int, sim_hours)) hours...")
progress = make_progress_callback(sim_hours, update_interval_hours=250)
@time sol = run_simulation(moon, sim_hours, T0, callback=progress)
println()  # newline after progress

println("\nSimulation complete!")
println("  Number of timesteps: $(length(sol.t))")
println("  Final time: $(sol.t[end]/3600) hours")

# Reshape final state to 2D for analysis
T_final = reshape(sol.u[end], moon.n_lat, moon.n_lon)
T_final_C = T_final .- 273.15

# Global statistics
global_mean = sum(T_final_C .* moon.cell_areas)
global_min = minimum(T_final_C)
global_max = maximum(T_final_C)

println("\n" * "="^60)
println("GLOBAL STATISTICS")
println("="^60)
println("  Mean temperature: $(round(global_mean, digits=1))°C")
println("  Min temperature:  $(round(global_min, digits=1))°C")
println("  Max temperature:  $(round(global_max, digits=1))°C")

# Find hottest and coldest spots
max_idx = argmax(T_final_C)
min_idx = argmin(T_final_C)
hot_lat = moon.latitudes[max_idx[1]]
hot_lon = moon.longitudes[max_idx[2]]
cold_lat = moon.latitudes[min_idx[1]]
cold_lon = moon.longitudes[min_idx[2]]

println("\n" * "="^60)
println("EXTREMES")
println("="^60)
println("  Hottest: $(round(global_max, digits=1))°C at lat=$(round(hot_lat, digits=1))°, lon=$(round(hot_lon, digits=1))°")
println("  Coldest: $(round(global_min, digits=1))°C at lat=$(round(cold_lat, digits=1))°, lon=$(round(cold_lon, digits=1))°")

# === Use unified output system ===
config = RunConfig("2d")
ctx = initialize_run(config)

save_results!(ctx, sol, moon, T0=T0)
generate_plots!(ctx, sol, moon)
generate_animations!(ctx, sol, moon, hours=448, frame_skip=4, fps=10)
finalize_run!(ctx)
