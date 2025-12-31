"""
Main script to run 2D Hot Moon simulation with two-layer atmosphere
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Statistics

include(joinpath(@__DIR__, "../src/HotMoon.jl"))
using .HotMoon

println("="^60)
println("Hot Moon 2D Climate Simulation (Two-Layer Atmosphere)")
println("="^60)

# Create moon with 45 lat × 90 lon grid (4° × 4° resolution)
println("\nCreating moon...")
moon = HotMoonBody(45, 90, seed=32, sea_level=0.0, scale=0.01, octaves=5)
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

# Initial upper layer state (uniform)
U0 = fill(U_INITIAL, moon.n_lat, moon.n_lon)
M_up0 = fill(M_UP_INITIAL, moon.n_lat, moon.n_lon)
println("Initial upper mass: $(U_INITIAL) (uniform)")
println("Initial upper moisture: $(M_UP_INITIAL) kg/m² (uniform)")

# Run simulation
sim_hours = 15000.0  # Long enough to see circulation develop
println("\nRunning two-layer atmosphere simulation for $(round(Int, sim_hours)) hours...")
progress = make_progress_callback(sim_hours, update_interval_hours=100)
@time sol = run_simulation_with_twolayer_atmosphere(moon, sim_hours, T0, M0, U0=U0, M_up0=M_up0, callback=progress)
println()  # newline after progress

println("\nSimulation complete!")
println("  Number of timesteps: $(length(sol.t))")
println("  Final time: $(sol.t[end]/3600) hours")

# Extract final state
n_cells = moon.n_lat * moon.n_lon
T_final = reshape(sol.u[end][1:n_cells], moon.n_lat, moon.n_lon)
M_final = reshape(sol.u[end][n_cells+1:2n_cells], moon.n_lat, moon.n_lon)
U_final = reshape(sol.u[end][2n_cells+1:3n_cells], moon.n_lat, moon.n_lon)
M_up_final = reshape(sol.u[end][3n_cells+1:4n_cells], moon.n_lat, moon.n_lon)
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
println("SURFACE MOISTURE STATISTICS")
println("="^60)
println("  Mean moisture: $(round(global_mean_M, digits=2)) kg/m²")
println("  Min moisture:  $(round(global_min_M, digits=2)) kg/m²")
println("  Max moisture:  $(round(global_max_M, digits=2)) kg/m²")

# Upper mass statistics
global_mean_U = sum(U_final .* moon.cell_areas)
global_min_U = minimum(U_final)
global_max_U = maximum(U_final)

println("\n" * "="^60)
println("UPPER MASS STATISTICS")
println("="^60)
println("  Mean U: $(round(global_mean_U, digits=3))")
println("  Min U:  $(round(global_min_U, digits=3))")
println("  Max U:  $(round(global_max_U, digits=3))")

# Upper moisture statistics
global_mean_M_up = sum(M_up_final .* moon.cell_areas)
global_min_M_up = minimum(M_up_final)
global_max_M_up = maximum(M_up_final)

println("\n" * "="^60)
println("UPPER MOISTURE STATISTICS")
println("="^60)
println("  Mean M_up: $(round(global_mean_M_up, digits=4)) kg/m²")
println("  Min M_up:  $(round(global_min_M_up, digits=4)) kg/m²")
println("  Max M_up:  $(round(global_max_M_up, digits=4)) kg/m²")

# Compute precipitation field at final time (including descent suppression)
precip_final = zeros(moon.n_lat, moon.n_lon)
for i in 1:moon.n_lat
    for j in 1:moon.n_lon
        lat = moon.latitudes[i]
        U_cell = max(U_FLOOR, U_final[i, j])
        descent = compute_total_descent_rate(U_cell, lat)
        drying_factor = compute_descent_saturation_multiplier(descent)
        M_sat_effective = compute_saturation_moisture_at_temperature(T_final[i, j]) * drying_factor
        if M_final[i, j] > M_sat_effective
            precip_final[i, j] = PRECIP_RATE * (M_final[i, j] - M_sat_effective)
        end
    end
end
precip_mm_hr = precip_final .* 3600

println("\n" * "="^60)
println("PRECIPITATION STATISTICS (with descent suppression)")
println("="^60)
println("  Mean precip rate: $(round(mean(precip_mm_hr), digits=3)) mm/hr")
println("  Max precip rate:  $(round(maximum(precip_mm_hr), digits=3)) mm/hr")
wet_cells = sum(precip_mm_hr .> 0.001)
wet_pct = 100 * wet_cells / n_cells
dry_cells = sum(precip_mm_hr .< 0.0001)
dry_pct = 100 * dry_cells / n_cells
println("  Cells with precipitation: $(round(wet_pct, digits=1))%")
println("  Cells effectively dry:    $(round(dry_pct, digits=1))%")

# Circulation diagnostics
println("\n" * "="^60)
println("CIRCULATION DIAGNOSTICS")
println("="^60)

# Zonal mean U by latitude
zonal_U = [mean(U_final[i, :]) for i in 1:moon.n_lat]
max_U_lat_idx = argmax(zonal_U)
min_U_lat_idx = argmin(zonal_U)
println("  Peak U latitude:   $(round(moon.latitudes[max_U_lat_idx], digits=1))° (U=$(round(zonal_U[max_U_lat_idx], digits=3)))")
println("  Minimum U latitude: $(round(moon.latitudes[min_U_lat_idx], digits=1))° (U=$(round(zonal_U[min_U_lat_idx], digits=3)))")

# Find descent zones (high U, low surface moisture)
println("\n  Looking for subtropical desert signatures...")
for lat_deg in [-30.0, -20.0, 20.0, 30.0]
    lat_idx = argmin(abs.(moon.latitudes .- lat_deg))
    avg_U = mean(U_final[lat_idx, :])
    avg_M = mean(M_final[lat_idx, :])
    avg_precip = mean(precip_mm_hr[lat_idx, :])
    println("    Lat $(round(Int, lat_deg))°: U=$(round(avg_U, digits=3)), M=$(round(avg_M, digits=2)) kg/m², precip=$(round(avg_precip, digits=3)) mm/hr")
end

# Mass conservation check
println("\n  Mass conservation check:")
println("    Initial total U: $(round(sum(U0 .* moon.cell_areas), digits=3))")
println("    Final total U:   $(round(sum(U_final .* moon.cell_areas), digits=3))")

# Print biome statistics
print_biome_statistics(T_final, M_final, moon)

# === Use unified output system ===
config = RunConfig("2d_twolayer")
ctx = initialize_run(config)

save_results!(ctx, sol, moon, T0=T0, M0=M0)
generate_plots!(ctx, sol, moon)
generate_animations!(ctx, sol, moon, hours=448, frame_skip=4, fps=10)
finalize_run!(ctx)
