"""
Plotting and visualization functions for 1D and 2D simulations
"""

"""
    plot_global_mean(sol, moon::MoonBody1D)

Plot global mean temperature over time.
"""
function plot_global_mean(sol, moon::MoonBody1D)
    n_points = length(sol.t)
    start_idx = max(1, n_points - 800)

    times_hr = sol.t[start_idx:end] ./ 3600
    global_temps = [sum(sol.u[i] .* moon.cell_areas) for i in start_idx:n_points]
    global_temps_C = global_temps .- 273.15

    p = plot(times_hr, global_temps_C,
        xlabel="Time (hours)",
        ylabel="Global Mean Temperature (°C)",
        title="Hot Moon - Global Mean Temperature",
        linewidth=2,
        legend=false,
        size=(800, 400))

    t_start_hr = times_hr[1]
    t_end_hr = times_hr[end]

    rotation_period_hr = ROTATION_PERIOD / 3600


    for cycle_start in 0:rotation_period_hr:(t_end_hr+rotation_period_hr)
        night_start = cycle_start + 0.25 * rotation_period_hr
        night_end = cycle_start + 0.75 * rotation_period_hr

        if night_end >= t_start_hr && night_start <= t_end_hr
            vspan!([night_start, night_end],
                color=:navy, alpha=0.08, label="")
        end

        day1_start = cycle_start
        day1_end = cycle_start + 0.25 * rotation_period_hr
        if day1_end >= t_start_hr && day1_start <= t_end_hr
            vspan!([day1_start, day1_end],
                color=:yellow, alpha=0.05, label="")
        end

        day2_start = cycle_start + 0.75 * rotation_period_hr
        day2_end = cycle_start + rotation_period_hr
        if day2_end >= t_start_hr && day2_start <= t_end_hr
            vspan!([day2_start, day2_end],
                color=:yellow, alpha=0.05, label="")
        end
    end

    # Add eclipse shading (on top of day/night)
    orbital_period_hr = ORBITAL_PERIOD / 3600
    eclipse_duration_hr = ECLIPSE_DURATION / 3600

    for t_hr in 0:orbital_period_hr:(t_end_hr+orbital_period_hr)
        eclipse_center = t_hr + orbital_period_hr / 2
        eclipse_start = eclipse_center - eclipse_duration_hr / 2
        eclipse_end = eclipse_center + eclipse_duration_hr / 2

        if eclipse_end >= t_start_hr && eclipse_start <= t_end_hr
            vspan!([eclipse_start, eclipse_end],
                color=:gray, alpha=0.3, label="")
        end
    end

    return p
end

"""
    plot_latitude_profile(sol, moon::MoonBody1D)

Plot temperature vs latitude (time-averaged).
"""
function plot_latitude_profile(sol, moon::MoonBody1D)
    # Average over last 20% of simulation (equilibrium)
    n_avg = div(length(sol.t), 5)
    temps_avg = mean([sol.u[i] for i in (length(sol.t)-n_avg):length(sol.t)])
    temps_avg_C = temps_avg .- 273.15

    plot(moon.latitudes, temps_avg_C,
        xlabel="Latitude (°)",
        ylabel="Temperature (°C)",
        title="Temperature Profile by Latitude",
        linewidth=2,
        marker=:circle,
        legend=false,
        size=(800, 400))

    # Add reference lines
    hline!([0], linestyle=:dash, color=:blue, label="Freezing", alpha=0.5)
end

"""
    plot_heatmap(sol, moon::MoonBody1D)

Create a heatmap of temperature evolution.
"""
function plot_heatmap(sol, moon::MoonBody1D)
    times_hr = sol.t ./ 3600

    # Convert to matrix: [latitude, time]
    temps_matrix = hcat(sol.u...)
    temps_matrix_C = temps_matrix .- 273.15

    heatmap(times_hr, moon.latitudes, temps_matrix_C,
        xlabel="Time (hours)",
        ylabel="Latitude (°)",
        title="Temperature Evolution",
        color=:thermal,
        size=(900, 400))
end

"""
    plot_latitude_mean_range(sol, moon::MoonBody1D)

Plot mean temperature at each latitude with min/max range bands.
Shows how much temperature varies (day/night, eclipses) at each latitude.
"""
function plot_latitude_mean_range(sol, moon::MoonBody1D)
    # Use last 20% of simulation (equilibrium)
    n_start = max(1, length(sol.t) - div(length(sol.t), 5))

    n_lat = moon.n_lat
    temps_mean = zeros(n_lat)
    temps_min = zeros(n_lat)
    temps_max = zeros(n_lat)

    for i in 1:n_lat
        band_temps = [sol.u[j][i] for j in n_start:length(sol.t)]
        temps_mean[i] = mean(band_temps)
        temps_min[i] = minimum(band_temps)
        temps_max[i] = maximum(band_temps)
    end

    # Convert to Celsius
    temps_mean_C = temps_mean .- 273.15
    temps_min_C = temps_min .- 273.15
    temps_max_C = temps_max .- 273.15

    p = plot(moon.latitudes, temps_mean_C,
        xlabel="Latitude (°)",
        ylabel="Temperature (°C)",
        title="Temperature by Latitude (Mean ± Range)",
        linewidth=3,
        color=:red,
        label="Mean",
        legend=:topright,
        size=(800, 500))

    # Add shaded range
    plot!(moon.latitudes, temps_min_C,
        fillrange=temps_max_C,
        fillalpha=0.3,
        fillcolor=:orange,
        linewidth=1,
        linestyle=:dash,
        color=:orange,
        label="Day/Night Range")

    # Reference lines
    hline!([0], linestyle=:dash, color=:blue, alpha=0.5, label="Freezing (0°C)")
    hline!([100], linestyle=:dot, color=:red, alpha=0.3, label="Boiling (100°C)")

    # Mark habitable zone (light green band from 0-40°C)
    hspan!([0, 40], color=:green, alpha=0.1, label="")

    return p
end

"""
    plot_latitude_timeseries(sol, moon::MoonBody1D)

Plot temperature over time for each latitude band.
Like global mean plot, but shows all latitudes as separate lines.
"""
function plot_latitude_timeseries(sol, moon::MoonBody1D)
    n_points = length(sol.t)
    start_idx = max(1, n_points - 1100)

    times_hr = sol.t[start_idx:end] ./ 3600

    # Create color gradient from warm (equator) to cool (pole)
    colors = cgrad(:thermal, moon.n_lat, categorical=true)

    p = plot(
        xlabel="Time (hours)",
        ylabel="Temperature (°C)",
        title="Temperature by Latitude Over Time",
        legend=:outerright,
        size=(1000, 500))

    # Plot each latitude band
    for i in 1:moon.n_lat
        temps_C = [sol.u[j][i] - 273.15 for j in start_idx:n_points]
        lat_label = "$(round(Int, moon.latitudes[i]))°"
        plot!(times_hr, temps_C,
            color=colors[i],
            linewidth=1.5,
            label=lat_label,
            alpha=0.8)
    end

    # Add eclipse shading
    t_start_hr = times_hr[1]
    t_end_hr = times_hr[end]
    orbital_period_hr = ORBITAL_PERIOD / 3600
    eclipse_duration_hr = ECLIPSE_DURATION / 3600

    for t_hr in 0:orbital_period_hr:(t_end_hr+orbital_period_hr)
        eclipse_center = t_hr + orbital_period_hr / 2
        eclipse_start = eclipse_center - eclipse_duration_hr / 2
        eclipse_end = eclipse_center + eclipse_duration_hr / 2

        if eclipse_end >= t_start_hr && eclipse_start <= t_end_hr
            vspan!([eclipse_start, eclipse_end],
                color=:gray, alpha=0.2, label="")
        end
    end

    # Reference line for freezing
    hline!([0], linestyle=:dash, color=:blue, alpha=0.5, label="")

    return p
end

"""
    plot_summary(sol, moon::MoonBody1D, filename="summary.png")

Create and save individual plots.
"""
function plot_summary(sol, moon::MoonBody1D, filename="summary.png")
    # Extract base name and directory
    base = replace(filename, r"\.(png|pdf|svg)$" => "")

    # Create and save each plot individually
    p1 = plot_global_mean(sol, moon)
    savefig(p1, "$(base)_global_mean.png")
    println("Saved: $(base)_global_mean.png")

    p2 = plot_latitude_profile(sol, moon)
    savefig(p2, "$(base)_latitude_profile.png")
    println("Saved: $(base)_latitude_profile.png")

    p3 = plot_heatmap(sol, moon)
    savefig(p3, "$(base)_heatmap.png")
    println("Saved: $(base)_heatmap.png")

    p4 = plot_latitude_mean_range(sol, moon)
    savefig(p4, "$(base)_latitude_range.png")
    println("Saved: $(base)_latitude_range.png")

    p5 = plot_latitude_timeseries(sol, moon)
    savefig(p5, "$(base)_latitude_timeseries.png")
    println("Saved: $(base)_latitude_timeseries.png")
end

# =============================================================================
# 2D Visualization Functions
# =============================================================================

"""
    add_eclipse_shading!(t_start_hr, t_end_hr)

Add vertical shading for eclipse periods to the current plot.
"""
function add_eclipse_shading!(t_start_hr, t_end_hr)
    orbital_period_hr = ORBITAL_PERIOD / 3600
    eclipse_duration_hr = ECLIPSE_DURATION / 3600

    for t_hr in 0:orbital_period_hr:(t_end_hr + orbital_period_hr)
        eclipse_center = t_hr + orbital_period_hr / 2
        eclipse_start = eclipse_center - eclipse_duration_hr / 2
        eclipse_end = eclipse_center + eclipse_duration_hr / 2

        if eclipse_end >= t_start_hr && eclipse_start <= t_end_hr
            vspan!([eclipse_start, eclipse_end], color=:gray, alpha=0.3, label="")
        end
    end
end

"""
    plot_global_mean_full(sol, moon::MoonBody2D)

Plot global mean temperature for the ENTIRE simulation.
Use this to check if the simulation has reached equilibrium.
"""
function plot_global_mean_full(sol, moon::MoonBody2D)
    times_hr = sol.t ./ 3600

    # Compute global mean at each timestep
    global_temps_C = Float64[]
    for i in 1:length(sol.t)
        T_2d = reshape(sol.u[i], moon.n_lat, moon.n_lon)
        push!(global_temps_C, sum(T_2d .* moon.cell_areas) - 273.15)
    end

    p = plot(times_hr, global_temps_C,
        xlabel="Time (hours)",
        ylabel="Global Mean Temperature (°C)",
        title="Global Mean Temperature - Full Simulation",
        linewidth=1.5,
        legend=false,
        color=:red,
        size=(1000, 400))

    return p
end

"""
    plot_global_mean_detail(sol, moon::MoonBody2D; hours=800)

Plot global mean temperature for the last N hours.
Shows eclipse impacts clearly with shading.
"""
function plot_global_mean_detail(sol, moon::MoonBody2D; hours=800)
    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    times_hr = sol.t[start_idx:end] ./ 3600

    # Compute global mean at each timestep
    global_temps_C = Float64[]
    for i in start_idx:length(sol.t)
        T_2d = reshape(sol.u[i], moon.n_lat, moon.n_lon)
        push!(global_temps_C, sum(T_2d .* moon.cell_areas) - 273.15)
    end

    p = plot(times_hr, global_temps_C,
        xlabel="Time (hours)",
        ylabel="Global Mean Temperature (°C)",
        title="Global Mean Temperature - Last $(hours) hours",
        linewidth=2,
        legend=false,
        color=:red,
        size=(1000, 400))

    add_eclipse_shading!(times_hr[1], times_hr[end])

    return p
end

"""
    plot_snapshot(sol, moon::MoonBody2D; time_idx=-1)

Single heatmap snapshot of temperature across the globe.
Uses percentile-based color limits to preserve detail in equatorial regions
despite very cold poles.
"""
function plot_snapshot(sol, moon::MoonBody2D; time_idx=-1)
    idx = time_idx < 0 ? length(sol.t) : time_idx
    T_2d = reshape(sol.u[idx], moon.n_lat, moon.n_lon) .- 273.15

    # Use 10th percentile as lower limit to not let poles dominate
    sorted_temps = sort(vec(T_2d))
    p10 = sorted_temps[max(1, div(length(sorted_temps), 10))]
    p100 = sorted_temps[end]

    time_hr = round(sol.t[idx] / 3600, digits=1)

    p = heatmap(moon.longitudes, moon.latitudes, T_2d,
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
        title="Temperature Snapshot (t = $(time_hr) hr)",
        color=:thermal,
        clims=(p10, p100),
        colorbar_title="°C",
        size=(1000, 500))

    return p
end

"""
    plot_hovmoeller_longitude(sol, moon::MoonBody2D; lat_idx=nothing, hours=280)

Time-Longitude (Hovmöller) diagram at a specific latitude (default: equator).
X-axis: time, Y-axis: longitude, color: temperature.
Shows the day/night terminator moving across longitudes over time.
Default 280 hours = 10 rotation periods (10 days).
"""
function plot_hovmoeller_longitude(sol, moon::MoonBody2D; lat_idx=nothing, hours=280)
    # Default to equator
    if lat_idx === nothing
        lat_idx = argmin(abs.(moon.latitudes))
    end

    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    n_times = length(sol.t) - start_idx + 1
    n_lon = moon.n_lon

    # Build matrix: rows = longitude, cols = time
    temps_matrix = zeros(n_lon, n_times)
    for (col, i) in enumerate(start_idx:length(sol.t))
        T_2d = reshape(sol.u[i], moon.n_lat, moon.n_lon)
        temps_matrix[:, col] = T_2d[lat_idx, :] .- 273.15
    end

    times_hr = sol.t[start_idx:end] ./ 3600
    lat_deg = round(moon.latitudes[lat_idx], digits=1)

    p = heatmap(times_hr, moon.longitudes, temps_matrix,
        xlabel="Time (hours)",
        ylabel="Longitude (°)",
        title="Time-Longitude at Latitude $(lat_deg)°",
        color=:thermal,
        colorbar_title="°C",
        size=(1000, 500))

    add_eclipse_shading!(times_hr[1], times_hr[end])

    return p
end

"""
    plot_hovmoeller_latitude(sol, moon::MoonBody2D; lon_idx=nothing, hours=280)

Time-Latitude (Hovmöller) diagram at a specific longitude (default: subsolar point, lon=0°).
X-axis: time, Y-axis: latitude, color: temperature.
Shows temperature evolution across latitudes over time.
Default 280 hours = 10 rotation periods (10 days).
"""
function plot_hovmoeller_latitude(sol, moon::MoonBody2D; lon_idx=nothing, hours=280)
    # Default to subsolar longitude (0°)
    if lon_idx === nothing
        lon_idx = argmin(abs.(moon.longitudes))
    end

    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    n_times = length(sol.t) - start_idx + 1
    n_lat = moon.n_lat

    # Build matrix: rows = latitude, cols = time
    temps_matrix = zeros(n_lat, n_times)
    for (col, i) in enumerate(start_idx:length(sol.t))
        T_2d = reshape(sol.u[i], moon.n_lat, moon.n_lon)
        temps_matrix[:, col] = T_2d[:, lon_idx] .- 273.15
    end

    times_hr = sol.t[start_idx:end] ./ 3600
    lon_deg = round(moon.longitudes[lon_idx], digits=1)

    p = heatmap(times_hr, moon.latitudes, temps_matrix,
        xlabel="Time (hours)",
        ylabel="Latitude (°)",
        title="Time-Latitude at Longitude $(lon_deg)°",
        color=:thermal,
        colorbar_title="°C",
        size=(1000, 500))

    add_eclipse_shading!(times_hr[1], times_hr[end])

    return p
end

"""
    plot_latitude_mean_range(sol, moon::MoonBody2D)

Plot temperature by latitude with mean and min/max range.
"""
function plot_latitude_mean_range(sol, moon::MoonBody2D)
    n_start = max(1, length(sol.t) - div(length(sol.t), 5))
    n_lat = moon.n_lat

    temps_mean = zeros(n_lat)
    temps_min = zeros(n_lat)
    temps_max = zeros(n_lat)

    for i in 1:n_lat
        all_temps = Float64[]
        for j in n_start:length(sol.t)
            T_2d = reshape(sol.u[j], moon.n_lat, moon.n_lon)
            append!(all_temps, T_2d[i, :])
        end
        temps_mean[i] = mean(all_temps) - 273.15
        temps_min[i] = minimum(all_temps) - 273.15
        temps_max[i] = maximum(all_temps) - 273.15
    end

    p = plot(moon.latitudes, temps_mean,
        xlabel="Latitude (°)",
        ylabel="Temperature (°C)",
        title="Temperature by Latitude (Mean + Range)",
        linewidth=3,
        color=:red,
        label="Mean",
        legend=:topright,
        size=(800, 500))

    plot!(moon.latitudes, temps_min,
        fillrange=temps_max,
        fillalpha=0.2,
        fillcolor=:orange,
        linewidth=1,
        linestyle=:dash,
        color=:orange,
        label="Min/Max Range")

    hline!([0], linestyle=:dash, color=:blue, alpha=0.5, label="")
    hspan!([0, 40], color=:green, alpha=0.1, label="")

    return p
end

"""
    plot_elevation_map(moon::MoonBody2D)

Create a heatmap of terrain elevation with coastline contour.
Ocean (elevation < 0) shown in blue tones, land in terrain colors.
"""
function plot_elevation_map(moon::MoonBody2D)
    p = heatmap(moon.longitudes, moon.latitudes, moon.elevation,
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
        title="Terrain Elevation Map",
        color=:terrain,
        colorbar_title="Elevation",
        size=(1000, 500))

    # Add coastline (elevation = 0)
    contour!(moon.longitudes, moon.latitudes, moon.elevation,
        levels=[0.0], color=:black, linewidth=1, label="")

    return p
end