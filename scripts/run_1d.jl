"""
Main script to run 1D Hot Moon simulation
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "../src/HotMoon.jl"))
using .HotMoon

println("="^60)
println("Hot Moon 1D Climate Simulation")
println("="^60)

# Create moon with 18 latitude bands
println("\nCreating moon...")
moon = HotMoonBody(18)
println(moon)

# Initial conditions: uniform 285K (12°C)
T0 = fill(285.0, moon.n_lat)
println("Initial temperature: $(T0[1] - 273.15)°C")

# Run simulation for 20000 hours
sim_hours = 20000.0
println("\nRunning simulation for $(round(Int, sim_hours)) hours...")
progress = make_progress_callback(sim_hours, update_interval_hours=500)
@time sol = run_simulation(moon, sim_hours, T0, callback=progress)
println()  # newline after progress

println("Simulation complete!")
println("  Number of timesteps: $(length(sol.t))")
println("  Final time: $(sol.t[end]/3600) hours")

# Calculate final statistics
final_temps = sol.u[end] .- 273.15
global_mean = sum(final_temps .* moon.cell_areas)
println("\nFinal statistics:")
println("  Global mean: $(round(global_mean, digits=1))°C")
println("  Equator: $(round(final_temps[1], digits=1))°C")
println("  Pole: $(round(final_temps[end], digits=1))°C")

# === Use unified output system ===
config = RunConfig("1d")
ctx = initialize_run(config)

save_results!(ctx, sol, moon, T0=T0)
generate_plots!(ctx, sol, moon)
finalize_run!(ctx)
