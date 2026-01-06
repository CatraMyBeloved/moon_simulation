"""
Animation functions for 2D climate visualizations.
"""

# =============================================================================
# Generic Animation Function
# =============================================================================

"""
    animate_field_evolution(sol, moon::MoonBody2D, ::Type{V}, filename; kwargs...) where V<:PlotVariable

Create an animated GIF of a variable over time.

# Arguments
- `sol`: ODE solution
- `moon`: MoonBody2D instance
- `V`: PlotVariable type (Temperature, Moisture, etc.)
- `filename`: Output filename for GIF

# Keywords
- `hours`: Number of hours to animate (default: 448 = 1 metacycle)
- `frame_skip`: Skip frames to reduce file size (default: 4)
- `fps`: Frames per second (default: 10)
"""
function animate_field_evolution(sol, moon::MoonBody2D, ::Type{V}, filename::String;
                                  hours=448, frame_skip=4, fps=10) where V<:PlotVariable
    start_idx = find_start_index_for_hours_before_end(sol, hours)
    frame_indices = start_idx:frame_skip:length(sol.t)

    println("    Creating $(variable_name(V)) animation: $(length(frame_indices)) frames...")

    # Compute global min/max for consistent color scale
    all_vals = Float64[]
    for idx in frame_indices
        field = extract_field_for_plotting(sol, moon, V, idx)
        # Filter out NaN/Inf values
        valid_vals = filter(isfinite, vec(field))
        append!(all_vals, valid_vals)
    end

    clims = if V == Temperature
        sorted_vals = sort(all_vals)
        n = length(sorted_vals)
        (sorted_vals[max(1, div(3 * n, 10))], sorted_vals[max(1, div(98 * n, 100))])
    else
        vmin, vmax = minimum(all_vals), maximum(all_vals)
        # Ensure minimum spread to avoid GKS color indexing errors
        if vmax - vmin < 1e-6
            vmax = vmin + 1.0
        end
        (vmin, vmax)
    end

    anim = @animate for (frame_num, idx) in enumerate(frame_indices)
        field = extract_field_for_plotting(sol, moon, V, idx)
        # Clamp field to clims range to avoid color indexing issues
        field = clamp.(field, clims[1], clims[2])
        time_hr = round(sol.t[idx] / 3600, digits=1)

        heatmap(moon.longitudes, moon.latitudes, field,
            xlabel="Longitude (°)",
            ylabel="Latitude (°)",
            title="$(variable_name(V)) (t = $(time_hr) hr)",
            color=heatmap_colorscheme(V),
            clims=clims,
            colorbar_title=colorbar_label(V),
            size=(900, 450))

        contour!(moon.longitudes, moon.latitudes, moon.elevation,
            levels=[0.0], color=:black, linewidth=1, label="")
    end

    gif(anim, filename, fps=fps)
    println("    Saved: $filename")
end

# Convenience wrappers
animate_temperature_evolution(sol, moon, filename; kwargs...) =
    animate_field_evolution(sol, moon, Temperature, filename; kwargs...)
animate_moisture_evolution(sol, moon, filename; kwargs...) =
    animate_field_evolution(sol, moon, Moisture, filename; kwargs...)
animate_precipitation_evolution(sol, moon, filename; kwargs...) =
    animate_field_evolution(sol, moon, Precipitation, filename; kwargs...)

# =============================================================================
# Combined Animation (T, M, Precip side by side)
# =============================================================================

"""
    animate_combined_climate(sol, moon::MoonBody2D, filename; kwargs...)

Create an animated GIF with Temperature, Moisture, and Precipitation side by side.

# Keywords
- `hours`: Number of hours to animate (default: 448 = 1 metacycle)
- `frame_skip`: Skip frames to reduce file size (default: 4)
- `fps`: Frames per second (default: 10)
"""
function animate_combined_climate(sol, moon::MoonBody2D, filename::String;
                                   hours=448, frame_skip=4, fps=10)
    start_idx = find_start_index_for_hours_before_end(sol, hours)
    frame_indices = start_idx:frame_skip:length(sol.t)

    println("    Creating combined climate animation: $(length(frame_indices)) frames...")

    # Compute global ranges for consistent color scales, filtering NaN/Inf
    T_vals, M_vals, P_vals = Float64[], Float64[], Float64[]
    for idx in frame_indices
        T_field = extract_field_for_plotting(sol, moon, Temperature, idx)
        M_field = extract_field_for_plotting(sol, moon, Moisture, idx)
        P_field = extract_field_for_plotting(sol, moon, Precipitation, idx)
        append!(T_vals, filter(isfinite, vec(T_field)))
        append!(M_vals, filter(isfinite, vec(M_field)))
        append!(P_vals, filter(isfinite, vec(P_field)))
    end

    # Temperature: percentile-based clims
    T_sorted = sort(T_vals)
    n_T = length(T_sorted)
    T_clims = (T_sorted[max(1, div(3 * n_T, 10))], T_sorted[max(1, div(98 * n_T, 100))])
    if T_clims[2] - T_clims[1] < 1e-6
        T_clims = (T_clims[1], T_clims[1] + 1.0)
    end

    # Moisture: min/max with minimum spread
    M_min, M_max = minimum(M_vals), maximum(M_vals)
    if M_max - M_min < 1e-6
        M_max = M_min + 1.0
    end
    M_clims = (M_min, M_max)

    # Precipitation: min/max with minimum spread
    P_min, P_max = minimum(P_vals), maximum(P_vals)
    if P_max - P_min < 1e-6
        P_max = P_min + 0.1
    end
    P_clims = (P_min, P_max)

    anim = @animate for (frame_num, idx) in enumerate(frame_indices)
        T_field = extract_field_for_plotting(sol, moon, Temperature, idx)
        M_field = extract_field_for_plotting(sol, moon, Moisture, idx)
        P_field = extract_field_for_plotting(sol, moon, Precipitation, idx)
        time_hr = round(sol.t[idx] / 3600, digits=1)

        # Clamp fields to avoid GKS color indexing errors
        T_field = clamp.(T_field, T_clims[1], T_clims[2])
        M_field = clamp.(M_field, M_clims[1], M_clims[2])
        P_field = clamp.(P_field, P_clims[1], P_clims[2])

        # Replace any remaining non-finite values with clims midpoint
        T_field = replace(T_field, NaN => (T_clims[1] + T_clims[2])/2, Inf => T_clims[2], -Inf => T_clims[1])
        M_field = replace(M_field, NaN => (M_clims[1] + M_clims[2])/2, Inf => M_clims[2], -Inf => M_clims[1])
        P_field = replace(P_field, NaN => (P_clims[1] + P_clims[2])/2, Inf => P_clims[2], -Inf => P_clims[1])

        p1 = heatmap(moon.longitudes, moon.latitudes, T_field,
            title="Temperature (°C) - t=$(time_hr)h",
            color=:thermal, clims=T_clims, colorbar=true)
        contour!(moon.longitudes, moon.latitudes, moon.elevation,
            levels=[0.0], color=:black, linewidth=0.5, label="")

        p2 = heatmap(moon.longitudes, moon.latitudes, M_field,
            title="Moisture (kg/m²)",
            color=:Blues, clims=M_clims, colorbar=true)
        contour!(moon.longitudes, moon.latitudes, moon.elevation,
            levels=[0.0], color=:black, linewidth=0.5, label="")

        p3 = heatmap(moon.longitudes, moon.latitudes, P_field,
            title="Precipitation (mm/hr)",
            color=:YlGnBu, clims=P_clims, colorbar=true)
        contour!(moon.longitudes, moon.latitudes, moon.elevation,
            levels=[0.0], color=:black, linewidth=0.5, label="")

        plot(p1, p2, p3, layout=(1, 3), size=(1500, 400))
    end

    gif(anim, filename, fps=fps)
    println("    Saved: $filename")
end

# =============================================================================
# Batch Animation
# =============================================================================

"""
    animate_all_variables(sol, moon::MoonBody2D, output_dir; kwargs...)

Generate all animation types for a moisture-coupled simulation.

# Keywords
- `hours`: Number of hours to animate (default: 448 = 1 metacycle)
- `frame_skip`: Skip frames to reduce file size (default: 4)
- `fps`: Frames per second (default: 10)
"""
function animate_all_variables(sol, moon::MoonBody2D, output_dir::String;
                                hours=448, frame_skip=4, fps=10)
    mkpath(output_dir)

    println("  Generating climate animation...")
    animate_combined_climate(sol, moon, joinpath(output_dir, "climate_animation.gif");
        hours=hours, frame_skip=frame_skip, fps=fps)

    println("  Generating temperature animation...")
    animate_field_evolution(sol, moon, Temperature,
        joinpath(output_dir, "temperature.gif"); hours=hours, frame_skip=frame_skip, fps=fps)

    println("  Generating moisture animation...")
    animate_field_evolution(sol, moon, Moisture,
        joinpath(output_dir, "moisture.gif"); hours=hours, frame_skip=frame_skip, fps=fps)

    println("  Generating precipitation animation...")
    animate_field_evolution(sol, moon, Precipitation,
        joinpath(output_dir, "precipitation.gif"); hours=hours, frame_skip=frame_skip, fps=fps)
end
