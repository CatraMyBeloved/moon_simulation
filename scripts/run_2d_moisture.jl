"""
Main script to run 2D Hot Moon simulation with moisture system
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Statistics

include(joinpath(@__DIR__, "../src/HotMoon.jl"))
using .HotMoon

println("="^60)
println("Hot Moon 2D Climate Simulation (with Moisture)")
println("="^60)

# Create moon with 90 lat × 180 lon grid (2° × 2° resolution)
println("\nCreating moon...")
moon = HotMoonBody(90, 180, seed=32, sea_level=0.0, scale=0.01, octaves=5)
println(moon)
println("  Total cells: $(moon.n_lat * moon.n_lon)")

# Calculate ocean coverage
ocean_cells = sum(moon.elevation .< 0)
ocean_pct = 100 * ocean_cells / (moon.n_lat * moon.n_lon)
println("  Ocean coverage: $(round(ocean_pct, digits=1))%")

# Initial conditions - use equilibrium estimate based on latitude/elevation
println("\nEstimating initial state from equilibrium conditions...")
T0, M0 = estimate_initial_temperature_and_moisture(moon)
println("Initial temperature range: $(round(minimum(T0) - 273.15, digits=1))°C to $(round(maximum(T0) - 273.15, digits=1))°C")
println("Initial moisture range: $(round(minimum(M0), digits=2)) to $(round(maximum(M0), digits=2)) kg/m²")

# Run simulation
sim_hours = 30000.0
println("\nRunning coupled T-M simulation for $(round(Int, sim_hours)) hours...")
progress = make_progress_callback(sim_hours, update_interval_hours=100)
@time sol = run_simulation_with_moisture(moon, sim_hours, T0, M0, callback=progress)
println()  # newline after progress

println("\nSimulation complete!")
println("  Number of timesteps: $(length(sol.t))")
println("  Final time: $(sol.t[end]/3600) hours")

# Extract final state
n_cells = moon.n_lat * moon.n_lon
T_final = reshape(sol.u[end][1:n_cells], moon.n_lat, moon.n_lon)
M_final = reshape(sol.u[end][n_cells+1:end], moon.n_lat, moon.n_lon)
T_final_C = T_final .- 273.15

# Global temperature statistics
global_mean_T = sum(T_final_C .* moon.cell_areas)
global_min_T = minimum(T_final_C)
global_max_T = maximum(T_final_C)

println("\n" * "="^60)
println("TEMPERATURE STATISTICS")
println("="^60)
println("  Mean temperature: $(round(global_mean_T, digits=1))°C")
println("  Min temperature:  $(round(global_min_T, digits=1))°C")
println("  Max temperature:  $(round(global_max_T, digits=1))°C")

# Global moisture statistics
global_mean_M = sum(M_final .* moon.cell_areas)
global_min_M = minimum(M_final)
global_max_M = maximum(M_final)

println("\n" * "="^60)
println("MOISTURE STATISTICS")
println("="^60)
println("  Mean moisture: $(round(global_mean_M, digits=2)) kg/m²")
println("  Min moisture:  $(round(global_min_M, digits=2)) kg/m²")
println("  Max moisture:  $(round(global_max_M, digits=2)) kg/m²")

# Compute precipitation field at final time
precip_final = zeros(moon.n_lat, moon.n_lon)
for i in 1:moon.n_lat
    for j in 1:moon.n_lon
        precip_final[i, j] = compute_orographic_precipitation_rate(M_final[i, j], T_final[i, j], moon.elevation[i, j])
    end
end
precip_mm_hr = precip_final .* 3600

println("\n" * "="^60)
println("PRECIPITATION STATISTICS")
println("="^60)
println("  Mean precip rate: $(round(mean(precip_mm_hr), digits=3)) mm/hr")
println("  Max precip rate:  $(round(maximum(precip_mm_hr), digits=3)) mm/hr")
wet_cells = sum(precip_mm_hr .> 0.001)
wet_pct = 100 * wet_cells / n_cells
println("  Cells with precipitation: $(round(wet_pct, digits=1))%")

# Find wettest and driest spots
max_M_idx = argmax(M_final)
min_M_idx = argmin(M_final)
wet_lat = moon.latitudes[max_M_idx[1]]
wet_lon = moon.longitudes[max_M_idx[2]]
dry_lat = moon.latitudes[min_M_idx[1]]
dry_lon = moon.longitudes[min_M_idx[2]]

println("\n" * "="^60)
println("MOISTURE EXTREMES")
println("="^60)
println("  Wettest: $(round(global_max_M, digits=2)) kg/m² at lat=$(round(wet_lat, digits=1))°, lon=$(round(wet_lon, digits=1))°")
println("  Driest:  $(round(global_min_M, digits=2)) kg/m² at lat=$(round(dry_lat, digits=1))°, lon=$(round(dry_lon, digits=1))°")

# Print biome statistics
print_biome_statistics(T_final, M_final, moon)

# === Use unified output system ===
config = RunConfig("2d_moisture")
ctx = initialize_run(config)

save_results!(ctx, sol, moon, T0=T0, M0=M0)
generate_plots!(ctx, sol, moon)
generate_animations!(ctx, sol, moon, hours=448, frame_skip=4, fps=10)
finalize_run!(ctx)
