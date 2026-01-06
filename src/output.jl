"""
Unified output management system for Hot Moon simulations.

Provides consistent run archiving, serialization, and visualization across all
simulation types (1D, 2D, 2D+moisture).

Note: Two-layer atmosphere output handling has been archived to src/archived/twolayer/

## Features
- Timestamped run directories
- Automatic copy to `latest/` folder
- Solution and moon serialization
- README.md generation with metadata
- Type-dispatched plot/animation generation

## Usage
```julia
config = RunConfig("2d_moisture")
ctx = initialize_run(config)
save_results!(ctx, sol, moon)
generate_plots!(ctx, sol, moon)
generate_animations!(ctx, sol, moon)
finalize_run!(ctx)
```
"""

using Dates
using Serialization
using Printf

# =============================================================================
# Configuration Types
# =============================================================================

"""
    RunConfig

Configuration for what to save and how.

# Fields
- `simulation_type`: One of "1d", "2d", "2d_moisture"
- `description`: Optional user description for the run
- `save_solution`: Whether to serialize the ODE solution
- `save_moon`: Whether to serialize the moon configuration
- `generate_plots`: Whether to create visualization PNGs
- `generate_animations`: Whether to create GIF animations (2D only)
"""
struct RunConfig
    simulation_type::String
    description::String
    save_solution::Bool
    save_moon::Bool
    generate_plots::Bool
    generate_animations::Bool
end

# Convenient constructors
RunConfig(sim_type::String) = RunConfig(sim_type, "", true, true, true, true)
RunConfig(sim_type::String, description::String) = RunConfig(sim_type, description, true, true, true, true)

function Base.show(io::IO, config::RunConfig)
    print(io, "RunConfig($(config.simulation_type))")
end

# =============================================================================
# Run Context (Internal State)
# =============================================================================

"""
    RunContext

Tracks state during a run for README generation and file management.
Created by `initialize_run()`, used by all output functions.
"""
mutable struct RunContext
    config::RunConfig
    run_dir::String           # Full path: output/runs/2024-12-31_193045_2d_moisture/
    plots_dir::String         # run_dir/plots/
    start_time::DateTime
    end_time::Union{DateTime,Nothing}

    # Tracked for README
    sim_hours::Float64
    n_timesteps::Int
    grid_size::Tuple
    ocean_pct::Float64
    initial_T_range::Tuple{Float64,Float64}
    final_T_range::Tuple{Float64,Float64}
    final_stats::Dict{String,Any}
    files_saved::Vector{String}

    # Terrain generation parameters (for reproducibility)
    terrain_seed::Union{Int,Nothing}
    terrain_sea_level::Union{Float64,Nothing}
    terrain_scale::Union{Float64,Nothing}
    terrain_octaves::Union{Int,Nothing}
end

function Base.show(io::IO, ctx::RunContext)
    print(io, "RunContext($(ctx.config.simulation_type) @ $(basename(ctx.run_dir)))")
end

# =============================================================================
# Run Initialization
# =============================================================================

"""
    initialize_run(config::RunConfig) -> RunContext

Create a timestamped run directory and return context for saving outputs.

Creates directory structure:
```
output/runs/YYYY-MM-DD_HHMMSS_<sim_type>/
└── plots/
```

# Example
```julia
ctx = initialize_run(RunConfig("2d_moisture"))
# ctx.run_dir = "output/runs/2024-12-31_193045_2d_moisture/"
```
"""
function initialize_run(config::RunConfig)
    # Create timestamped directory name
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    run_name = "$(timestamp)_$(config.simulation_type)"
    run_dir = joinpath("output", "runs", run_name)
    plots_dir = joinpath(run_dir, "plots")

    # Create directories
    mkpath(run_dir)
    mkpath(plots_dir)

    println("\n" * "="^60)
    println("INITIALIZING RUN: $run_name")
    println("="^60)
    println("  Directory: $run_dir")

    # Create context with default values
    ctx = RunContext(
        config,
        run_dir,
        plots_dir,
        now(),
        nothing,
        0.0,            # sim_hours
        0,              # n_timesteps
        (0, 0),         # grid_size
        0.0,            # ocean_pct
        (0.0, 0.0),     # initial_T_range
        (0.0, 0.0),     # final_T_range
        Dict{String,Any}(),
        String[],
        nothing,        # terrain_seed
        nothing,        # terrain_sea_level
        nothing,        # terrain_scale
        nothing         # terrain_octaves
    )

    return ctx
end

# =============================================================================
# Save Results
# =============================================================================

"""
    save_results!(ctx::RunContext, sol, moon; T0=nothing, M0=nothing, terrain_params=nothing, kwargs...)

Save solution and moon to the run directory, and record metadata.

# Arguments
- `ctx`: Run context from `initialize_run()`
- `sol`: ODE solution object
- `moon`: Moon body (MoonBody1D or MoonBody2D)
- `T0`: Initial temperature field (optional, for metadata)
- `M0`: Initial moisture field (optional, for metadata)
- `terrain_params`: NamedTuple with (seed, sea_level, scale, octaves) for reproducibility

Saves:
- `solution.jls`: Serialized ODE solution
- `moon.jls`: Serialized moon configuration
"""
function save_results!(ctx::RunContext, sol, moon; T0=nothing, M0=nothing, terrain_params=nothing, kwargs...)
    println("\n  Saving results...")

    # Record simulation metadata
    ctx.sim_hours = sol.t[end] / 3600
    ctx.n_timesteps = length(sol.t)

    # Grid info based on moon type
    if moon isa MoonBody1D
        ctx.grid_size = (moon.n_lat,)
        ctx.ocean_pct = 0.0
    else
        ctx.grid_size = (moon.n_lat, moon.n_lon)
        ctx.ocean_pct = 100 * sum(moon.elevation .< 0) / (moon.n_lat * moon.n_lon)
    end

    # Store terrain generation parameters if provided
    if terrain_params !== nothing
        ctx.terrain_seed = get(terrain_params, :seed, nothing)
        ctx.terrain_sea_level = get(terrain_params, :sea_level, nothing)
        ctx.terrain_scale = get(terrain_params, :scale, nothing)
        ctx.terrain_octaves = get(terrain_params, :octaves, nothing)
    end

    # Temperature ranges
    if T0 !== nothing
        T0_flat = T0 isa AbstractMatrix ? vec(T0) : T0
        ctx.initial_T_range = (minimum(T0_flat) - 273.15, maximum(T0_flat) - 273.15)
    end

    # Final temperature - extract from solution based on type
    T_final = extract_temperature_field(sol.u[end], moon)
    ctx.final_T_range = (minimum(T_final) - 273.15, maximum(T_final) - 273.15)

    # Save solution
    if ctx.config.save_solution
        sol_path = joinpath(ctx.run_dir, "solution.jls")
        serialize(sol_path, sol)
        push!(ctx.files_saved, "solution.jls")
        sol_size = filesize(sol_path) / 1024 / 1024
        println("    - solution.jls ($(round(sol_size, digits=1)) MB)")
    end

    # Save moon
    if ctx.config.save_moon
        moon_path = joinpath(ctx.run_dir, "moon.jls")
        serialize(moon_path, moon)
        push!(ctx.files_saved, "moon.jls")
        moon_size = filesize(moon_path) / 1024 / 1024
        println("    - moon.jls ($(round(moon_size, digits=1)) MB)")
    end

    # Store any additional stats
    for (k, v) in kwargs
        ctx.final_stats[string(k)] = v
    end
end

"""
Extract temperature field from solution state, handling different simulation types.
"""
function extract_temperature_field(u, moon::MoonBody1D)
    return u  # 1D: state is just temperature
end

function extract_temperature_field(u, moon::MoonBody2D)
    n_cells = moon.n_lat * moon.n_lon

    # Check state size to determine type
    if length(u) == n_cells
        # Temperature only
        return reshape(u, moon.n_lat, moon.n_lon)
    elseif length(u) == 2 * n_cells
        # Temperature + Moisture
        return reshape(u[1:n_cells], moon.n_lat, moon.n_lon)
    else
        error("Unknown state size: $(length(u))")
    end
end

# =============================================================================
# Plot Generation
# =============================================================================

"""
    generate_plots!(ctx::RunContext, sol, moon)

Generate all appropriate plots for the simulation type.
Dispatches to the correct plot set based on `ctx.config.simulation_type`.
"""
function generate_plots!(ctx::RunContext, sol, moon)
    if !ctx.config.generate_plots
        println("  Skipping plots (disabled in config)")
        return
    end

    println("\n" * "="^60)
    println("GENERATING PLOTS")
    println("="^60)

    sim_type = ctx.config.simulation_type

    if sim_type == "1d"
        _generate_plots_1d!(ctx, sol, moon)
    elseif sim_type == "2d"
        _generate_plots_2d!(ctx, sol, moon)
    elseif sim_type == "2d_moisture"
        _generate_plots_2d_moisture!(ctx, sol, moon)
    else
        @warn "Unknown simulation type: $sim_type"
    end
end

function _generate_plots_1d!(ctx::RunContext, sol, moon)
    println("  Creating 1D summary plots...")
    plot_1d_summary_to_files(sol, moon, joinpath(ctx.plots_dir, "summary.png"))
    push!(ctx.files_saved, "plots/summary.png")
end

function _generate_plots_2d!(ctx::RunContext, sol, moon)
    plots_dir = ctx.plots_dir

    # Temperature plots
    println("  Temperature plots...")
    _save_plot(plot_global_mean_timeseries(sol, moon, Temperature; hours=nothing),
               plots_dir, "temperature_global_mean_full.png", ctx)
    _save_plot(plot_global_mean_timeseries(sol, moon, Temperature; hours=800),
               plots_dir, "temperature_global_mean_detail.png", ctx)
    _save_plot(plot_field_snapshot(sol, moon, Temperature),
               plots_dir, "temperature.png", ctx)
    _save_plot(plot_longitude_time_hovmoeller(sol, moon, Temperature),
               plots_dir, "temperature_hovmoeller_longitude.png", ctx)
    _save_plot(plot_latitude_time_hovmoeller(sol, moon, Temperature),
               plots_dir, "temperature_hovmoeller_latitude.png", ctx)
    _save_plot(plot_latitude_mean_range(sol, moon, Temperature),
               plots_dir, "temperature_latitude_range.png", ctx)

    # Terrain
    println("  Terrain map...")
    _save_plot(plot_elevation_map(moon), plots_dir, "terrain.png", ctx)
end

function _generate_plots_2d_moisture!(ctx::RunContext, sol, moon)
    plots_dir = ctx.plots_dir

    # Temperature plots
    println("  Temperature plots...")
    _save_plot(plot_global_mean_timeseries(sol, moon, Temperature; hours=nothing),
               plots_dir, "temperature_global_mean_full.png", ctx)
    _save_plot(plot_global_mean_timeseries(sol, moon, Temperature; hours=800),
               plots_dir, "temperature_global_mean_detail.png", ctx)
    _save_plot(plot_field_snapshot(sol, moon, Temperature),
               plots_dir, "temperature.png", ctx)
    _save_plot(plot_longitude_time_hovmoeller(sol, moon, Temperature),
               plots_dir, "temperature_hovmoeller_longitude.png", ctx)
    _save_plot(plot_latitude_time_hovmoeller(sol, moon, Temperature),
               plots_dir, "temperature_hovmoeller_latitude.png", ctx)
    _save_plot(plot_latitude_mean_range(sol, moon, Temperature),
               plots_dir, "temperature_latitude_range.png", ctx)

    # Moisture plots
    println("  Moisture plots...")
    _save_plot(plot_global_mean_timeseries(sol, moon, Moisture; hours=nothing),
               plots_dir, "moisture_global_mean_full.png", ctx)
    _save_plot(plot_global_mean_timeseries(sol, moon, Moisture; hours=800),
               plots_dir, "moisture_global_mean_detail.png", ctx)
    _save_plot(plot_field_snapshot(sol, moon, Moisture),
               plots_dir, "moisture.png", ctx)
    _save_plot(plot_longitude_time_hovmoeller(sol, moon, Moisture),
               plots_dir, "moisture_hovmoeller_longitude.png", ctx)
    _save_plot(plot_latitude_time_hovmoeller(sol, moon, Moisture),
               plots_dir, "moisture_hovmoeller_latitude.png", ctx)
    _save_plot(plot_latitude_mean_range(sol, moon, Moisture),
               plots_dir, "moisture_latitude_range.png", ctx)

    # Precipitation plots
    println("  Precipitation plots...")
    _save_plot(plot_global_mean_timeseries(sol, moon, Precipitation; hours=nothing),
               plots_dir, "precipitation_global_mean_full.png", ctx)
    _save_plot(plot_global_mean_timeseries(sol, moon, Precipitation; hours=800),
               plots_dir, "precipitation_global_mean_detail.png", ctx)
    _save_plot(plot_field_snapshot(sol, moon, Precipitation),
               plots_dir, "precipitation.png", ctx)
    _save_plot(plot_longitude_time_hovmoeller(sol, moon, Precipitation),
               plots_dir, "precipitation_hovmoeller_longitude.png", ctx)
    _save_plot(plot_latitude_time_hovmoeller(sol, moon, Precipitation),
               plots_dir, "precipitation_hovmoeller_latitude.png", ctx)

    # Terrain and biome
    println("  Terrain and biome maps...")
    _save_plot(plot_elevation_map(moon), plots_dir, "terrain.png", ctx)
    _save_plot(plot_field_snapshot(sol, moon, Biome), plots_dir, "biome.png", ctx)
    _save_plot(plot_biome_with_legend(sol, moon), plots_dir, "biome_with_legend.png", ctx)
end

function _save_plot(p, dir, filename, ctx)
    path = joinpath(dir, filename)
    savefig(p, path)
    push!(ctx.files_saved, "plots/$filename")
end

# =============================================================================
# Animation Generation
# =============================================================================

"""
    generate_animations!(ctx::RunContext, sol, moon; hours=448, fps=10, frame_skip=4)

Generate animations for 2D simulations. No-op for 1D.

# Arguments
- `hours`: Duration to animate (default: 1 metacycle = 448 hours)
- `fps`: Frames per second in output GIF
- `frame_skip`: Skip frames to reduce file size
"""
function generate_animations!(ctx::RunContext, sol, moon; hours::Real=448, fps::Int=10, frame_skip::Int=4)
    if !ctx.config.generate_animations
        println("  Skipping animations (disabled in config)")
        return
    end

    if moon isa MoonBody1D
        println("  Skipping animations (1D simulation)")
        return
    end

    println("\n" * "="^60)
    println("GENERATING ANIMATIONS")
    println("="^60)

    sim_type = ctx.config.simulation_type
    plots_dir = ctx.plots_dir

    println("  Creating animated GIFs ($(hours) hours)...")

    if sim_type == "2d"
        animate_field_evolution(sol, moon, Temperature, joinpath(plots_dir, "temperature.gif"),
                               hours=hours, frame_skip=frame_skip, fps=fps)
        push!(ctx.files_saved, "plots/temperature.gif")

    elseif sim_type == "2d_moisture"
        animate_all_variables(sol, moon, plots_dir, hours=hours, frame_skip=frame_skip, fps=fps)
        push!(ctx.files_saved, "plots/temperature.gif")
        push!(ctx.files_saved, "plots/moisture.gif")
        push!(ctx.files_saved, "plots/precipitation.gif")
        push!(ctx.files_saved, "plots/biome.gif")
    end
end

# =============================================================================
# Finalize Run
# =============================================================================

"""
    finalize_run!(ctx::RunContext)

Complete the run: write README.md, copy files to `latest/`, print summary.

Must be called after all other output operations.
"""
function finalize_run!(ctx::RunContext)
    ctx.end_time = now()

    # Write README
    _write_readme(ctx)

    # Copy to latest/
    _copy_to_latest(ctx)

    # Print summary
    duration_sec = Dates.value(ctx.end_time - ctx.start_time) / 1000

    println("\n" * "="^60)
    println("RUN COMPLETE")
    println("="^60)
    println("  Run directory: $(ctx.run_dir)")
    println("  Duration: $(round(duration_sec, digits=1)) seconds")
    println("  Files saved: $(length(ctx.files_saved))")
    println("  Latest copy: output/latest/")
    println("="^60)
end

function _write_readme(ctx::RunContext)
    readme_path = joinpath(ctx.run_dir, "README.md")

    duration_sec = Dates.value(ctx.end_time - ctx.start_time) / 1000
    grid_str = length(ctx.grid_size) == 1 ? "$(ctx.grid_size[1]) latitude bands" :
               "$(ctx.grid_size[1]) lat × $(ctx.grid_size[2]) lon"

    open(readme_path, "w") do io
        println(io, "# Simulation Run: $(uppercase(ctx.config.simulation_type))")
        println(io)
        println(io, "**Date:** $(Dates.format(ctx.start_time, "yyyy-mm-dd HH:MM:SS"))")
        println(io, "**Duration:** $(round(duration_sec, digits=1)) seconds")
        println(io, "**Simulation type:** $(ctx.config.simulation_type)")
        println(io)

        println(io, "## Grid")
        println(io, "- Resolution: $grid_str")
        if length(ctx.grid_size) == 2
            println(io, "- Total cells: $(ctx.grid_size[1] * ctx.grid_size[2])")
            println(io, "- Ocean coverage: $(round(ctx.ocean_pct, digits=1))%")
        end
        println(io)

        # Terrain generation parameters (for reproducibility)
        if ctx.terrain_seed !== nothing
            println(io, "## Terrain Generation")
            println(io, "- seed: $(ctx.terrain_seed)")
            if ctx.terrain_sea_level !== nothing
                println(io, "- sea_level: $(ctx.terrain_sea_level)")
            end
            if ctx.terrain_scale !== nothing
                println(io, "- scale: $(ctx.terrain_scale)")
            end
            if ctx.terrain_octaves !== nothing
                println(io, "- octaves: $(ctx.terrain_octaves)")
            end
            println(io)
        end

        println(io, "## Simulation")
        println(io, "- Simulated time: $(round(ctx.sim_hours, digits=1)) hours")
        println(io, "- Timesteps saved: $(ctx.n_timesteps)")
        if ctx.initial_T_range != (0.0, 0.0)
            println(io, "- Initial T range: $(round(ctx.initial_T_range[1], digits=1))°C to $(round(ctx.initial_T_range[2], digits=1))°C")
        end
        println(io, "- Final T range: $(round(ctx.final_T_range[1], digits=1))°C to $(round(ctx.final_T_range[2], digits=1))°C")
        println(io)

        if !isempty(ctx.final_stats)
            println(io, "## Statistics")
            for (k, v) in ctx.final_stats
                println(io, "- $k: $v")
            end
            println(io)
        end

        # Count files
        n_plots = count(f -> endswith(f, ".png"), ctx.files_saved)
        n_anims = count(f -> endswith(f, ".gif"), ctx.files_saved)
        n_data = count(f -> endswith(f, ".jls"), ctx.files_saved)

        println(io, "## Files")
        println(io, "- $n_data data files (.jls)")
        println(io, "- $n_plots plots (.png)")
        println(io, "- $n_anims animations (.gif)")
        println(io)

        if !isempty(ctx.config.description)
            println(io, "## Description")
            println(io, ctx.config.description)
            println(io)
        end

        println(io, "---")
        println(io, "*Generated by HotMoon.jl*")
    end

    push!(ctx.files_saved, "README.md")
    println("  Wrote README.md")
end

function _copy_to_latest(ctx::RunContext)
    latest_dir = joinpath("output", "latest")

    # Remove old latest if exists
    if isdir(latest_dir)
        rm(latest_dir, recursive=true, force=true)
    end

    # Small delay to let Windows release file handles from GIF creation
    sleep(0.5)

    # Copy entire run directory to latest with retry logic for Windows
    try
        cp(ctx.run_dir, latest_dir, force=true)
        println("  Copied to output/latest/")
    catch e
        if Sys.iswindows()
            # Windows file locking - retry after longer delay
            sleep(1.0)
            try
                cp(ctx.run_dir, latest_dir, force=true)
                println("  Copied to output/latest/ (after retry)")
            catch e2
                @warn "Could not copy to latest/ (Windows file locking)" exception=e2
                println("  Warning: Could not copy to output/latest/ - files may be locked")
            end
        else
            rethrow(e)
        end
    end
end

# =============================================================================
# Loading and Listing Runs
# =============================================================================

"""
    RunInfo

Summary information about a saved run.
"""
struct RunInfo
    path::String
    simulation_type::String
    timestamp::DateTime
    has_solution::Bool
    has_moon::Bool
end

function Base.show(io::IO, info::RunInfo)
    print(io, "RunInfo($(info.simulation_type), $(Dates.format(info.timestamp, "yyyy-mm-dd HH:MM")))")
end

"""
    load_run(run_path::String) -> (sol, moon, metadata)

Reload a previous run from its directory.

# Arguments
- `run_path`: Path to run directory (e.g., "output/runs/2024-12-31_193045_2d_moisture")

# Returns
- `sol`: Deserialized ODE solution (or `nothing` if not saved)
- `moon`: Deserialized moon configuration (or `nothing` if not saved)
- `metadata`: Dict with information parsed from README.md

# Example
```julia
sol, moon, meta = load_run("output/latest")
```
"""
function load_run(run_path::String)
    # Expand relative paths
    if !isabspath(run_path)
        run_path = abspath(run_path)
    end

    if !isdir(run_path)
        error("Run directory not found: $run_path")
    end

    sol = nothing
    moon = nothing
    metadata = Dict{String,Any}()

    # Load solution if exists
    sol_path = joinpath(run_path, "solution.jls")
    if isfile(sol_path)
        println("Loading solution.jls...")
        sol = deserialize(sol_path)
    end

    # Load moon if exists
    moon_path = joinpath(run_path, "moon.jls")
    if isfile(moon_path)
        println("Loading moon.jls...")
        moon = deserialize(moon_path)
    end

    # Parse README for metadata
    readme_path = joinpath(run_path, "README.md")
    if isfile(readme_path)
        metadata["readme"] = read(readme_path, String)
        # Extract simulation type from first line
        lines = split(metadata["readme"], "\n")
        if length(lines) > 0
            m = match(r"Run: (.+)$", lines[1])
            if m !== nothing
                metadata["simulation_type"] = lowercase(strip(m.captures[1]))
            end
        end
    end

    metadata["path"] = run_path

    return (sol, moon, metadata)
end

"""
    list_runs() -> Vector{RunInfo}

List all saved runs with their metadata.

# Example
```julia
runs = list_runs()
for r in runs
    println("\$(r.timestamp): \$(r.simulation_type)")
end
```
"""
function list_runs()
    runs_dir = joinpath("output", "runs")
    if !isdir(runs_dir)
        return RunInfo[]
    end

    runs = RunInfo[]

    for name in readdir(runs_dir)
        run_path = joinpath(runs_dir, name)
        if !isdir(run_path)
            continue
        end

        # Parse timestamp and type from directory name
        # Format: YYYY-MM-DD_HHMMSS_<type>
        m = match(r"^(\d{4}-\d{2}-\d{2})_(\d{6})_(.+)$", name)
        if m === nothing
            continue
        end

        date_str = m.captures[1]
        time_str = m.captures[2]
        sim_type = m.captures[3]

        timestamp = try
            DateTime("$(date_str)T$(time_str[1:2]):$(time_str[3:4]):$(time_str[5:6])")
        catch
            continue
        end

        has_sol = isfile(joinpath(run_path, "solution.jls"))
        has_moon = isfile(joinpath(run_path, "moon.jls"))

        push!(runs, RunInfo(run_path, sim_type, timestamp, has_sol, has_moon))
    end

    # Sort by timestamp, newest first
    sort!(runs, by=r -> r.timestamp, rev=true)

    return runs
end

"""
    print_runs()

Print a formatted list of all saved runs.
"""
function print_runs()
    runs = list_runs()

    if isempty(runs)
        println("No saved runs found in output/runs/")
        return
    end

    println("="^60)
    println("SAVED RUNS")
    println("="^60)

    for (i, r) in enumerate(runs)
        time_str = Dates.format(r.timestamp, "yyyy-mm-dd HH:MM:SS")
        data_str = r.has_solution ? "✓ sol" : "✗ sol"
        data_str *= r.has_moon ? ", ✓ moon" : ", ✗ moon"

        println("  [$i] $time_str  $(rpad(r.simulation_type, 12))  $data_str")
    end

    println("="^60)
end
