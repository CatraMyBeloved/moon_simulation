"""
Main script to run 2D Hot Moon simulation with moisture system
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Statistics
using Plots

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

# Initial conditions
T0 = fill(285.0, moon.n_lat, moon.n_lon)  # 12°C uniform
M0 = fill(2.0, moon.n_lat, moon.n_lon)    # 2 kg/m² moisture everywhere
println("Initial temperature: $(T0[1,1] - 273.15)°C (uniform)")
println("Initial moisture: $(M0[1,1]) kg/m² (uniform)")

# Run simulation
sim_hours = 5000.0  # Shorter for testing moisture dynamics
println("\nRunning coupled T-M simulation for $(round(Int, sim_hours)) hours...")
progress = make_progress_callback(sim_hours, update_interval_hours=250)
@time sol = run_simulation_moisture(moon, sim_hours, T0, M0, callback=progress)
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
        precip_final[i, j] = get_precipitation(M_final[i, j], T_final[i, j], moon.elevation[i, j])
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

# Generate visualizations
println("\n" * "="^60)
println("GENERATING PLOTS")
println("="^60)

mkpath("output/plots_2d_moisture")

println("  Creating temperature snapshot...")
p_temp = plot_temperature_snapshot_moisture(sol, moon)
savefig(p_temp, "output/plots_2d_moisture/temperature.png")

println("  Creating moisture snapshot...")
p_moist = plot_moisture_snapshot(sol, moon)
savefig(p_moist, "output/plots_2d_moisture/moisture.png")

println("  Creating precipitation map...")
p_precip = plot_precipitation_snapshot(sol, moon)
savefig(p_precip, "output/plots_2d_moisture/precipitation.png")

println("  Creating terrain map...")
p_terrain = plot_elevation_map(moon)
savefig(p_terrain, "output/plots_2d_moisture/terrain.png")

println("\n" * "="^60)
println("Done! Plots saved to output/plots_2d_moisture/")
println("="^60)
