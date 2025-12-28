"""
Load and analyze saved simulation results
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Serialization
using Statistics

include(joinpath(@__DIR__, "../src/HotMoon.jl"))
using .HotMoon

# Get the run directory from command line or use most recent
if length(ARGS) > 0
    run_dir = ARGS[1]
else
    # Find most recent run
    runs_dir = "output/runs"
    if isdir(runs_dir)
        run_folders = filter(d -> isdir(joinpath(runs_dir, d)), readdir(runs_dir))
        if !isempty(run_folders)
            run_dir = joinpath(runs_dir, last(sort(run_folders)))
            println("Loading most recent run: $run_dir")
        else
            error("No saved runs found in $runs_dir")
        end
    else
        error("No runs directory found")
    end
end

# Load the data
println("Loading simulation data from $run_dir...")
sol = deserialize("$run_dir/solution.jls")
moon = deserialize("$run_dir/moon.jls")

println("Data loaded successfully!")
println("  Simulation time: $(sol.t[end]/3600) hours")
println("  Timesteps: $(length(sol.t))")
println("  Latitude bands: $(moon.n_lat)")

# Calculate statistics based on model type
if moon isa MoonBody1D
    final_temps = sol.u[end] .- 273.15
    global_mean = sum(final_temps .* moon.cell_areas)
    println("\nFinal statistics:")
    println("  Global mean: $(round(global_mean, digits=1))°C")
    println("  Equator: $(round(final_temps[1], digits=1))°C")
    println("  Pole: $(round(final_temps[end], digits=1))°C")

    # Regenerate plots
    println("\nRegenerating plots...")
    plot_summary(sol, moon, "$run_dir/reanalyzed.png")
else
    # MoonBody2D
    T_2d = reshape(sol.u[end], moon.n_lat, moon.n_lon) .- 273.15
    global_mean = sum(T_2d .* moon.cell_areas)
    equator_idx = argmin(abs.(moon.latitudes))

    println("\nFinal statistics:")
    println("  Global mean: $(round(global_mean, digits=1))°C")
    println("  Min temperature: $(round(minimum(T_2d), digits=1))°C")
    println("  Max temperature: $(round(maximum(T_2d), digits=1))°C")
    println("  Equator mean: $(round(mean(T_2d[equator_idx, :]), digits=1))°C")

    # Regenerate plots
    println("\nRegenerating plots...")
    p1 = plot_global_mean_full(sol, moon)
    savefig(p1, "$run_dir/reanalyzed_global_mean.png")
    p2 = plot_snapshot(sol, moon)
    savefig(p2, "$run_dir/reanalyzed_snapshot.png")
    println("Saved: $(run_dir)/reanalyzed_global_mean.png")
    println("Saved: $(run_dir)/reanalyzed_snapshot.png")
end

println("\nDone! New plots saved to $run_dir")
println("\nYou can now interactively analyze the data:")
println("  - sol.t contains time points (seconds)")
println("  - sol.u contains temperature arrays (Kelvin) at each time")
println("  - moon contains the spatial grid information")
