"""
Main script to run 1D Hot Moon simulation
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using Serialization

include(joinpath(@__DIR__, "../src/HotMoon.jl"))
using .HotMoon

println("="^60)
println("Hot Moon 1D Climate Simulation")
println("="^60)

mkpath(joinpath(@__DIR__, "../output/plots"))
mkpath(joinpath(@__DIR__, "../output/animations"))
mkpath(joinpath(@__DIR__, "../output/data"))

# Create moon with 18 latitude bands
println("\nCreating moon...")
moon = HotMoonBody(18)
println(moon)

# Initial conditions: uniform 285K (12°C)
T0 = fill(285.0, moon.n_lat)
println("Initial temperature: $(T0[1] - 273.15)°C")

# Run simulation for 2000 hours
println("\nRunning simulation for 2000 hours...")
@time sol = run_simulation(moon, 20000.0, T0)

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

# Create plots
println("\nCreating plots...")

# Save to timestamped directory
timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
run_dir = "output/runs/run_$timestamp"
mkpath(run_dir)

# Save plots with timestamp
plot_summary(sol, moon, "$run_dir/hotmoon_1d.png")

plot_summary(sol, moon, "output/plots/hotmoon_1d.png")

println("\n" * "="^60)
println("Done! Check $run_dir or output/plots/hotmoon_1d.png")
println("="^60)