"""
Main script to run 2D Hot Moon simulation
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Statistics
using Plots

include(joinpath(@__DIR__, "../src/HotMoon.jl"))
using .HotMoon

println("="^60)
println("Hot Moon 2D Climate Simulation")
println("="^60)

# Create moon with 36 lat × 72 lon grid (2.5° × 5° resolution)
println("\nCreating moon...")
moon = HotMoonBody(90, 180, seed=32, sea_level=0.05, scale=0.014, octaves=5)
println(moon)
println("  Total cells: $(moon.n_lat * moon.n_lon)")

# Initial conditions: uniform 285K (12°C)
T0 = fill(285.0, moon.n_lat, moon.n_lon)
println("Initial temperature: $(T0[1,1] - 273.15)°C (uniform)")

# Run simulation
sim_hours = 10000.0
println("\nRunning simulation for $(sim_hours) hours...")
println("  (This may take a while...)")
@time sol = run_simulation(moon, sim_hours, T0)

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
# Latitude profile (average across longitudes)
println("\n" * "="^60)
println("LATITUDE PROFILE (averaged across longitudes)")
println("="^60)
for i in 1:moon.n_lat
    lat = moon.latitudes[i]
    T_mean = mean(T_final_C[i, :])
    T_min = minimum(T_final_C[i, :])
    T_max = maximum(T_final_C[i, :])
    println("  $(lpad(round(Int, lat), 2))°: mean=$(lpad(round(T_mean, digits=1), 6))°C, " *
            "range=[$(lpad(round(T_min, digits=1), 6)), $(lpad(round(T_max, digits=1), 6))]°C")
end

# Longitude slice at equator
println("\n" * "="^60)
println("EQUATOR LONGITUDE PROFILE (lat ≈ $(round(moon.latitudes[1], digits=1))°)")
println("="^60)
equator_temps = T_final_C[1, :]
for j in 1:4:moon.n_lon  # Every 4th longitude to keep output manageable
    lon = moon.longitudes[j]
    T = equator_temps[j]
    bar = repeat("█", max(0, round(Int, (T + 50) / 5)))  # Simple ASCII bar
    println("  $(lpad(round(Int, lon), 3))°: $(lpad(round(T, digits=1), 6))°C  $bar")
end

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

# Generate visualizations
println("\n" * "="^60)
println("GENERATING PLOTS")
println("="^60)

mkpath("output/plots_2d")

println("  Creating global mean (full simulation)...")
p1 = plot_global_mean_full(sol, moon)
savefig(p1, "output/plots_2d/global_mean_full.png")

println("  Creating global mean (last 800h detail)...")
p2 = plot_global_mean_detail(sol, moon, hours=800)
savefig(p2, "output/plots_2d/global_mean_detail.png")

println("  Creating temperature snapshot...")
p3 = plot_snapshot(sol, moon)
savefig(p3, "output/plots_2d/snapshot.png")

println("  Creating Hovmöller diagram...")
p4 = plot_hovmoeller(sol, moon, hours=400)
savefig(p4, "output/plots_2d/hovmoeller.png")

println("  Creating latitude profile...")
p5 = plot_latitude_mean_range(sol, moon)
savefig(p5, "output/plots_2d/latitude_range.png")

println("\n" * "="^60)
println("Done! Plots saved to output/plots_2d/")
println("="^60)
