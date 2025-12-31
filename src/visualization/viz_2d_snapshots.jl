"""
2D snapshot plots: heatmap visualizations at a specific time.
"""

# =============================================================================
# Field Snapshots
# =============================================================================

"""
    plot_field_snapshot(sol, moon::MoonBody2D, variable; time_idx=-1, show_coastlines=true)

Single heatmap snapshot at a specific time.
time_idx=-1 shows final state.
"""
function plot_field_snapshot(sol, moon::MoonBody2D, ::Type{V};
                              time_idx=-1, show_coastlines=true) where V<:PlotVariable
    idx = time_idx < 0 ? length(sol.t) : time_idx
    field = extract_field_for_plotting(sol, moon, V, idx)

    time_hr = round(sol.t[idx] / 3600, digits=1)

    # Compute color limits
    clims = compute_snapshot_color_limits(field, V)

    p = heatmap(moon.longitudes, moon.latitudes, field,
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
        title="$(variable_name(V)) Snapshot (t = $(time_hr) hr)",
        color=heatmap_colorscheme(V),
        clims=clims,
        colorbar_title=colorbar_label(V),
        size=(1000, 500))

    if show_coastlines
        add_coastline_contour_to_plot!(moon)
    end

    return p
end

function compute_snapshot_color_limits(field, ::Type{Temperature})
    sorted_vals = sort(vec(field))
    n = length(sorted_vals)
    p30 = sorted_vals[max(1, div(3 * n, 10))]
    p98 = sorted_vals[max(1, div(98 * n, 100))]
    return (p30, p98)
end

function compute_snapshot_color_limits(field, ::Type{V}) where V<:PlotVariable
    return (minimum(field), maximum(field))
end

# Convenience wrappers
plot_temperature_snapshot(sol, moon; kwargs...) = plot_field_snapshot(sol, moon, Temperature; kwargs...)
plot_moisture_snapshot(sol, moon; kwargs...) = plot_field_snapshot(sol, moon, Moisture; kwargs...)
plot_precipitation_snapshot(sol, moon; kwargs...) = plot_field_snapshot(sol, moon, Precipitation; kwargs...)
plot_upper_mass_snapshot(sol, moon; kwargs...) = plot_field_snapshot(sol, moon, UpperMass; kwargs...)
plot_upper_moisture_snapshot(sol, moon; kwargs...) = plot_field_snapshot(sol, moon, UpperMoisture; kwargs...)

# =============================================================================
# Terrain Map
# =============================================================================

"""
    plot_elevation_map(moon::MoonBody2D)

Create a heatmap of terrain elevation with coastline contour.
"""
function plot_elevation_map(moon::MoonBody2D)
    p = heatmap(moon.longitudes, moon.latitudes, moon.elevation,
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
        title="Terrain Elevation Map",
        color=:terrain,
        colorbar_title="Elevation",
        size=(1000, 500))

    add_coastline_contour_to_plot!(moon)

    return p
end

# =============================================================================
# Biome Visualization
# =============================================================================

"""
    compute_biome_map(T_2d, M_2d, moon::MoonBody2D) -> Matrix{Int}

Compute biome classification for each cell from temperature, moisture, and elevation.
Returns matrix of biome IDs (0-11).
"""
function compute_biome_map(T_2d, M_2d, moon::MoonBody2D)
    biome_map = zeros(Int, moon.n_lat, moon.n_lon)
    for i in 1:moon.n_lat
        for j in 1:moon.n_lon
            biome_map[i, j] = classify_cell_biome(T_2d[i, j], M_2d[i, j], moon.elevation[i, j])
        end
    end
    return biome_map
end

"""
    plot_field_snapshot(sol, moon::MoonBody2D, ::Type{Biome}; time_idx=-1)

Create a biome classification map at a specific time.
Uses categorical colors with a legend showing biome names.
"""
function plot_field_snapshot(sol, moon::MoonBody2D, ::Type{Biome}; time_idx=-1, show_coastlines=true)
    idx = time_idx < 0 ? length(sol.t) : time_idx
    T_2d = extract_temperature_field_celsius(sol, moon, idx) .+ 273.15
    M_2d = extract_moisture_field(sol, moon, idx)
    biome_map = compute_biome_map(T_2d, M_2d, moon)

    time_hr = round(sol.t[idx] / 3600, digits=1)

    biome_cmap = cgrad([RGB(c...) for c in BIOME_COLORS_RGB], NUM_BIOMES, categorical=true)

    p = heatmap(moon.longitudes, moon.latitudes, biome_map,
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
        title="Biome Classification (t = $(time_hr) hr)",
        color=biome_cmap,
        clims=(-0.5, NUM_BIOMES - 0.5),
        colorbar=true,
        size=(1000, 500))

    if show_coastlines
        add_coastline_contour_to_plot!(moon)
    end

    return p
end

"""
    print_biome_statistics(T_2d, M_2d, moon::MoonBody2D)

Print statistics about biome coverage.
Shows percentage of surface area covered by each biome type.
"""
function print_biome_statistics(T_2d, M_2d, moon::MoonBody2D)
    biome_map = compute_biome_map(T_2d, M_2d, moon)

    biome_areas = zeros(NUM_BIOMES)
    for i in 1:moon.n_lat
        for j in 1:moon.n_lon
            biome_id = biome_map[i, j]
            biome_areas[biome_id + 1] += moon.cell_areas[i, j]
        end
    end
    total_area = sum(moon.cell_areas)

    println("\nBiome Coverage (area-weighted):")
    println("-" ^ 40)
    for id in 0:(NUM_BIOMES-1)
        pct = 100 * biome_areas[id + 1] / total_area
        if pct > 0.1
            println(@sprintf("  %-20s: %5.1f%%", BIOME_NAMES[id + 1], pct))
        end
    end
    println("-" ^ 40)
end

"""
    plot_biome_with_legend(sol, moon::MoonBody2D; time_idx=-1)

Create a biome map with a separate legend panel showing all biome types.
"""
function plot_biome_with_legend(sol, moon::MoonBody2D; time_idx=-1)
    idx = time_idx < 0 ? length(sol.t) : time_idx
    T_2d = extract_temperature_field_celsius(sol, moon, idx) .+ 273.15
    M_2d = extract_moisture_field(sol, moon, idx)
    biome_map = compute_biome_map(T_2d, M_2d, moon)

    time_hr = round(sol.t[idx] / 3600, digits=1)

    biome_cmap = cgrad([RGB(c...) for c in BIOME_COLORS_RGB], NUM_BIOMES, categorical=true)

    p1 = heatmap(moon.longitudes, moon.latitudes, biome_map,
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
        title="Biome Classification (t = $(time_hr) hr)",
        color=biome_cmap,
        clims=(-0.5, NUM_BIOMES - 0.5),
        colorbar=false,
        size=(800, 400))

    add_coastline_contour_to_plot!(moon)

    # Create legend as scatter plot with colored markers
    legend_x = fill(0.5, NUM_BIOMES)
    legend_y = collect(NUM_BIOMES:-1:1)

    p2 = scatter(legend_x, legend_y,
        marker=:square,
        markersize=15,
        markercolor=[RGB(c...) for c in BIOME_COLORS_RGB],
        markerstrokewidth=0,
        xlim=(0, 5),
        ylim=(0, NUM_BIOMES + 1),
        axis=false,
        grid=false,
        legend=false,
        size=(200, 400))

    for (i, name) in enumerate(BIOME_NAMES)
        annotate!(p2, 1.0, NUM_BIOMES - i + 1, text(name, :left, 8))
    end

    plot(p1, p2, layout=@layout([a{0.8w} b{0.2w}]), size=(1200, 450))
end
