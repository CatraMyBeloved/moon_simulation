"""
Quick script to visualize terrain generation
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Plots

include(joinpath(@__DIR__, "../src/HotMoon.jl"))
using .HotMoon

println("Generating terrain...")
moon = HotMoonBody(90, 180, seed=32, sea_level=0.05, scale=0.01, octaves=5)

println("Creating elevation map...")
p = plot_elevation_map(moon)

mkpath("output")
savefig(p, "output/terrain_map.png")
println("Saved to output/terrain_map.png")

# Print some stats
elev = moon.elevation
n_total = moon.n_lat * moon.n_lon
n_ocean = sum(elev .< 0)
n_land = sum(elev .>= 0)
println("\nTerrain Statistics:")
println("  Ocean cells: $n_ocean ($(round(100*n_ocean/n_total, digits=1))%)")
println("  Land cells:  $n_land ($(round(100*n_land/n_total, digits=1))%)")
println("  Elevation range: $(round(minimum(elev), digits=2)) to $(round(maximum(elev), digits=2))")
