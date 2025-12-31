"""
2D time series plots: global means and Hovmöller diagrams.
Uses generic functions to avoid duplication across variable types.
"""

# =============================================================================
# Global Mean Plots
# =============================================================================

"""
    plot_global_mean_timeseries(sol, moon::MoonBody2D, variable; hours=nothing, show_eclipses=true)

Plot global mean for any variable over time.
If hours=nothing, shows full simulation; otherwise shows last N hours.
"""
function plot_global_mean_timeseries(sol, moon::MoonBody2D, ::Type{V};
                                      hours=nothing, show_eclipses=true) where V<:PlotVariable
    if hours === nothing
        start_idx = 1
        title_suffix = "Full Simulation"
    else
        start_idx = find_start_index_for_hours_before_end(sol, hours)
        title_suffix = "Last $(hours) hours"
    end

    times_hr = convert_time_to_hours(sol, start_idx)

    global_vals = Float64[]
    for i in start_idx:length(sol.t)
        field = extract_field_for_plotting(sol, moon, V, i)
        push!(global_vals, sum(field .* moon.cell_areas))
    end

    p = plot(times_hr, global_vals,
        xlabel="Time (hours)",
        ylabel="Global Mean $(axis_label(V))",
        title="Global Mean $(variable_name(V)) - $(title_suffix)",
        linewidth=hours === nothing ? 1.5 : 2,
        legend=false,
        color=plot_line_color(V),
        size=(1000, 400))

    if show_eclipses && length(times_hr) > 1
        add_eclipse_shading_to_plot!(times_hr[1], times_hr[end])
    end

    return p
end

# Convenience wrappers
plot_temperature_global_mean(sol, moon; kwargs...) = plot_global_mean_timeseries(sol, moon, Temperature; kwargs...)
plot_moisture_global_mean(sol, moon; kwargs...) = plot_global_mean_timeseries(sol, moon, Moisture; kwargs...)
plot_precipitation_global_mean(sol, moon; kwargs...) = plot_global_mean_timeseries(sol, moon, Precipitation; kwargs...)
plot_upper_mass_global_mean(sol, moon; kwargs...) = plot_global_mean_timeseries(sol, moon, UpperMass; kwargs...)
plot_upper_moisture_global_mean(sol, moon; kwargs...) = plot_global_mean_timeseries(sol, moon, UpperMoisture; kwargs...)

# =============================================================================
# Hovmöller Diagrams
# =============================================================================

"""
    plot_longitude_time_hovmoeller(sol, moon::MoonBody2D, variable; latitude_index=nothing, hours=280)

Time-Longitude Hovmöller diagram at a specific latitude.
Default latitude is equator (latitude closest to 0°).
"""
function plot_longitude_time_hovmoeller(sol, moon::MoonBody2D, ::Type{V};
                                         latitude_index=nothing, hours=280) where V<:PlotVariable
    lat_idx = latitude_index === nothing ? argmin(abs.(moon.latitudes)) : latitude_index

    start_idx = find_start_index_for_hours_before_end(sol, hours)
    n_times = length(sol.t) - start_idx + 1

    data_matrix = zeros(moon.n_lon, n_times)
    for (col, i) in enumerate(start_idx:length(sol.t))
        field = extract_field_for_plotting(sol, moon, V, i)
        data_matrix[:, col] = field[lat_idx, :]
    end

    times_hr = convert_time_to_hours(sol, start_idx)
    lat_deg = round(moon.latitudes[lat_idx], digits=1)

    p = heatmap(times_hr, moon.longitudes, data_matrix,
        xlabel="Time (hours)",
        ylabel="Longitude (°)",
        title="Time-Longitude $(variable_name(V)) at Lat $(lat_deg)°",
        color=heatmap_colorscheme(V),
        colorbar_title=colorbar_label(V),
        size=(1000, 500))

    add_eclipse_shading_to_plot!(times_hr[1], times_hr[end])
    return p
end

"""
    plot_latitude_time_hovmoeller(sol, moon::MoonBody2D, variable; longitude_index=nothing, hours=280)

Time-Latitude Hovmöller diagram at a specific longitude.
Default longitude is the subsolar point (longitude closest to 0°).
"""
function plot_latitude_time_hovmoeller(sol, moon::MoonBody2D, ::Type{V};
                                        longitude_index=nothing, hours=280) where V<:PlotVariable
    lon_idx = longitude_index === nothing ? argmin(abs.(moon.longitudes)) : longitude_index

    start_idx = find_start_index_for_hours_before_end(sol, hours)
    n_times = length(sol.t) - start_idx + 1

    data_matrix = zeros(moon.n_lat, n_times)
    for (col, i) in enumerate(start_idx:length(sol.t))
        field = extract_field_for_plotting(sol, moon, V, i)
        data_matrix[:, col] = field[:, lon_idx]
    end

    times_hr = convert_time_to_hours(sol, start_idx)
    lon_deg = round(moon.longitudes[lon_idx], digits=1)

    p = heatmap(times_hr, moon.latitudes, data_matrix,
        xlabel="Time (hours)",
        ylabel="Latitude (°)",
        title="Time-Latitude $(variable_name(V)) at Lon $(lon_deg)°",
        color=heatmap_colorscheme(V),
        colorbar_title=colorbar_label(V),
        size=(1000, 500))

    add_eclipse_shading_to_plot!(times_hr[1], times_hr[end])
    return p
end

# =============================================================================
# Latitude Mean + Range Plots
# =============================================================================

"""
    plot_latitude_mean_range(sol, moon::MoonBody2D, variable)

Plot variable by latitude with mean and min/max range.
"""
function plot_latitude_mean_range(sol, moon::MoonBody2D, ::Type{V}) where V<:PlotVariable
    n_start = max(1, length(sol.t) - div(length(sol.t), 5))
    n_lat = moon.n_lat

    vals_mean = zeros(n_lat)
    vals_min = zeros(n_lat)
    vals_max = zeros(n_lat)

    for i in 1:n_lat
        all_vals = Float64[]
        for j in n_start:length(sol.t)
            field = extract_field_for_plotting(sol, moon, V, j)
            append!(all_vals, field[i, :])
        end
        vals_mean[i] = mean(all_vals)
        vals_min[i] = minimum(all_vals)
        vals_max[i] = maximum(all_vals)
    end

    p = plot(moon.latitudes, vals_mean,
        xlabel="Latitude (°)",
        ylabel=axis_label(V),
        title="$(variable_name(V)) by Latitude (Mean + Range)",
        linewidth=3,
        color=plot_line_color(V),
        label="Mean",
        legend=:topright,
        size=(800, 500))

    plot!(moon.latitudes, vals_min,
        fillrange=vals_max,
        fillalpha=0.2,
        fillcolor=plot_line_color(V),
        linewidth=1,
        linestyle=:dash,
        color=plot_line_color(V),
        alpha=0.5,
        label="Min/Max Range")

    # Add reference lines for temperature
    if V == Temperature
        hline!([0], linestyle=:dash, color=:blue, alpha=0.5, label="")
        hspan!([0, 40], color=:green, alpha=0.1, label="")
    end

    return p
end
