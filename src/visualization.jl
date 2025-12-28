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

# -----------------------------------------------------------------------------
# Plot Variable Types for Multiple Dispatch
# -----------------------------------------------------------------------------

"""
Abstract type for plot variables. Concrete subtypes enable dispatch-based
selection of what to plot from simulation results.
"""
abstract type PlotVariable end

"Temperature variable (in Kelvin internally, displayed as °C)"
struct Temperature <: PlotVariable end

"Moisture/water content variable (kg/m²)"
struct Moisture <: PlotVariable end

"Precipitation rate (computed from T and M, displayed as mm/hr)"
struct Precipitation <: PlotVariable end

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

"""
    _is_coupled_solution(sol, moon::MoonBody2D) -> Bool

Detect if solution is from a coupled T+M simulation based on state vector size.
"""
function _is_coupled_solution(sol, moon::MoonBody2D)
    n_cells = moon.n_lat * moon.n_lon
    return length(sol.u[1]) == 2 * n_cells
end

"""
    _extract_temperature(sol, moon::MoonBody2D, idx::Int) -> Matrix

Extract temperature field from solution at time index idx.
Auto-detects coupled vs temperature-only solutions.
"""
function _extract_temperature(sol, moon::MoonBody2D, idx::Int)
    n_cells = moon.n_lat * moon.n_lon
    if _is_coupled_solution(sol, moon)
        return reshape(sol.u[idx][1:n_cells], moon.n_lat, moon.n_lon)
    else
        return reshape(sol.u[idx], moon.n_lat, moon.n_lon)
    end
end

"""
    _extract_moisture(sol, moon::MoonBody2D, idx::Int) -> Matrix

Extract moisture field from coupled solution at time index idx.
"""
function _extract_moisture(sol, moon::MoonBody2D, idx::Int)
    n_cells = moon.n_lat * moon.n_lon
    return reshape(sol.u[idx][n_cells+1:end], moon.n_lat, moon.n_lon)
end

"""
    _compute_precipitation(T_2d, M_2d, moon::MoonBody2D) -> Matrix

Compute precipitation field from temperature and moisture.
Returns precipitation in kg/m²/s.
"""
function _compute_precipitation(T_2d, M_2d, moon::MoonBody2D)
    precip = zeros(moon.n_lat, moon.n_lon)
    for i in 1:moon.n_lat
        for j in 1:moon.n_lon
            precip[i, j] = get_precipitation(M_2d[i, j], T_2d[i, j], moon.elevation[i, j])
        end
    end
    return precip
end

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

# -----------------------------------------------------------------------------
# Plot Configuration by Variable Type
# -----------------------------------------------------------------------------

_plot_color(::Type{Temperature}) = :red
_plot_color(::Type{Moisture}) = :blue
_plot_color(::Type{Precipitation}) = :teal

_heatmap_color(::Type{Temperature}) = :thermal
_heatmap_color(::Type{Moisture}) = :Blues
_heatmap_color(::Type{Precipitation}) = :YlGnBu

_ylabel(::Type{Temperature}) = "Temperature (°C)"
_ylabel(::Type{Moisture}) = "Moisture (kg/m²)"
_ylabel(::Type{Precipitation}) = "Precipitation (mm/hr)"

_colorbar_title(::Type{Temperature}) = "°C"
_colorbar_title(::Type{Moisture}) = "kg/m²"
_colorbar_title(::Type{Precipitation}) = "mm/hr"

_variable_name(::Type{Temperature}) = "Temperature"
_variable_name(::Type{Moisture}) = "Moisture"
_variable_name(::Type{Precipitation}) = "Precipitation"

# -----------------------------------------------------------------------------
# Global Mean Plots
# -----------------------------------------------------------------------------

"""
    plot_global_mean_full(sol, moon::MoonBody2D, [variable])

Plot global mean for the ENTIRE simulation.
- `variable`: Temperature (default), Moisture, or Precipitation
"""
function plot_global_mean_full(sol, moon::MoonBody2D, ::Type{Temperature}=Temperature)
    times_hr = sol.t ./ 3600

    global_vals = Float64[]
    for i in 1:length(sol.t)
        T_2d = _extract_temperature(sol, moon, i)
        push!(global_vals, sum(T_2d .* moon.cell_areas) - 273.15)
    end

    p = plot(times_hr, global_vals,
        xlabel="Time (hours)",
        ylabel="Global Mean Temperature (°C)",
        title="Global Mean Temperature - Full Simulation",
        linewidth=1.5,
        legend=false,
        color=_plot_color(Temperature),
        size=(1000, 400))

    return p
end

function plot_global_mean_full(sol, moon::MoonBody2D, ::Type{Moisture})
    times_hr = sol.t ./ 3600

    global_vals = Float64[]
    for i in 1:length(sol.t)
        M_2d = _extract_moisture(sol, moon, i)
        push!(global_vals, sum(M_2d .* moon.cell_areas))
    end

    p = plot(times_hr, global_vals,
        xlabel="Time (hours)",
        ylabel="Global Mean Moisture (kg/m²)",
        title="Global Mean Moisture - Full Simulation",
        linewidth=1.5,
        legend=false,
        color=_plot_color(Moisture),
        size=(1000, 400))

    return p
end

function plot_global_mean_full(sol, moon::MoonBody2D, ::Type{Precipitation})
    times_hr = sol.t ./ 3600

    global_vals = Float64[]
    for i in 1:length(sol.t)
        T_2d = _extract_temperature(sol, moon, i)
        M_2d = _extract_moisture(sol, moon, i)
        precip = _compute_precipitation(T_2d, M_2d, moon)
        push!(global_vals, sum(precip .* moon.cell_areas) * 3600)  # to mm/hr
    end

    p = plot(times_hr, global_vals,
        xlabel="Time (hours)",
        ylabel="Global Mean Precipitation (mm/hr)",
        title="Global Mean Precipitation - Full Simulation",
        linewidth=1.5,
        legend=false,
        color=_plot_color(Precipitation),
        size=(1000, 400))

    return p
end

"""
    plot_global_mean_detail(sol, moon::MoonBody2D, [variable]; hours=800)

Plot global mean for the last N hours with eclipse shading.
- `variable`: Temperature (default), Moisture, or Precipitation
"""
function plot_global_mean_detail(sol, moon::MoonBody2D, ::Type{Temperature}=Temperature; hours=800)
    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    times_hr = sol.t[start_idx:end] ./ 3600

    global_vals = Float64[]
    for i in start_idx:length(sol.t)
        T_2d = _extract_temperature(sol, moon, i)
        push!(global_vals, sum(T_2d .* moon.cell_areas) - 273.15)
    end

    p = plot(times_hr, global_vals,
        xlabel="Time (hours)",
        ylabel="Global Mean Temperature (°C)",
        title="Global Mean Temperature - Last $(hours) hours",
        linewidth=2,
        legend=false,
        color=_plot_color(Temperature),
        size=(1000, 400))

    add_eclipse_shading!(times_hr[1], times_hr[end])
    return p
end

function plot_global_mean_detail(sol, moon::MoonBody2D, ::Type{Moisture}; hours=800)
    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    times_hr = sol.t[start_idx:end] ./ 3600

    global_vals = Float64[]
    for i in start_idx:length(sol.t)
        M_2d = _extract_moisture(sol, moon, i)
        push!(global_vals, sum(M_2d .* moon.cell_areas))
    end

    p = plot(times_hr, global_vals,
        xlabel="Time (hours)",
        ylabel="Global Mean Moisture (kg/m²)",
        title="Global Mean Moisture - Last $(hours) hours",
        linewidth=2,
        legend=false,
        color=_plot_color(Moisture),
        size=(1000, 400))

    add_eclipse_shading!(times_hr[1], times_hr[end])
    return p
end

function plot_global_mean_detail(sol, moon::MoonBody2D, ::Type{Precipitation}; hours=800)
    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    times_hr = sol.t[start_idx:end] ./ 3600

    global_vals = Float64[]
    for i in start_idx:length(sol.t)
        T_2d = _extract_temperature(sol, moon, i)
        M_2d = _extract_moisture(sol, moon, i)
        precip = _compute_precipitation(T_2d, M_2d, moon)
        push!(global_vals, sum(precip .* moon.cell_areas) * 3600)
    end

    p = plot(times_hr, global_vals,
        xlabel="Time (hours)",
        ylabel="Global Mean Precipitation (mm/hr)",
        title="Global Mean Precipitation - Last $(hours) hours",
        linewidth=2,
        legend=false,
        color=_plot_color(Precipitation),
        size=(1000, 400))

    add_eclipse_shading!(times_hr[1], times_hr[end])
    return p
end

# -----------------------------------------------------------------------------
# Snapshot Plots
# -----------------------------------------------------------------------------

"""
    plot_snapshot(sol, moon::MoonBody2D, [variable]; time_idx=-1)

Single heatmap snapshot at a specific time.
- `variable`: Temperature (default), Moisture, or Precipitation
- `time_idx`: Time index (-1 for final state)
"""
function plot_snapshot(sol, moon::MoonBody2D, ::Type{Temperature}=Temperature; time_idx=-1)
    idx = time_idx < 0 ? length(sol.t) : time_idx
    T_2d = _extract_temperature(sol, moon, idx) .- 273.15

    sorted_vals = sort(vec(T_2d))
    n = length(sorted_vals)
    p30 = sorted_vals[max(1, div(3 * n, 10))]
    p98 = sorted_vals[max(1, div(98 * n, 100))]

    time_hr = round(sol.t[idx] / 3600, digits=1)

    p = heatmap(moon.longitudes, moon.latitudes, T_2d,
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
        title="Temperature Snapshot (t = $(time_hr) hr)",
        color=:inferno,
        clims=(p30, p98),
        colorbar_title="°C",
        size=(1000, 500))

    return p
end

function plot_snapshot(sol, moon::MoonBody2D, ::Type{Moisture}; time_idx=-1)
    idx = time_idx < 0 ? length(sol.t) : time_idx
    M_2d = _extract_moisture(sol, moon, idx)

    time_hr = round(sol.t[idx] / 3600, digits=1)

    p = heatmap(moon.longitudes, moon.latitudes, M_2d,
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
        title="Moisture (t = $(time_hr) hr)",
        color=:Blues,
        colorbar_title="kg/m²",
        size=(1000, 500))

    contour!(moon.longitudes, moon.latitudes, moon.elevation,
        levels=[0.0], color=:black, linewidth=1, label="")

    return p
end

function plot_snapshot(sol, moon::MoonBody2D, ::Type{Precipitation}; time_idx=-1)
    idx = time_idx < 0 ? length(sol.t) : time_idx
    T_2d = _extract_temperature(sol, moon, idx)
    M_2d = _extract_moisture(sol, moon, idx)
    precip_mm_hr = _compute_precipitation(T_2d, M_2d, moon) .* 3600

    time_hr = round(sol.t[idx] / 3600, digits=1)

    p = heatmap(moon.longitudes, moon.latitudes, precip_mm_hr,
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
        title="Precipitation Rate (t = $(time_hr) hr)",
        color=:YlGnBu,
        colorbar_title="mm/hr",
        size=(1000, 500))

    contour!(moon.longitudes, moon.latitudes, moon.elevation,
        levels=[0.0], color=:black, linewidth=1, label="")

    return p
end

# -----------------------------------------------------------------------------
# Hovmöller Diagrams
# -----------------------------------------------------------------------------

"""
    plot_hovmoeller_longitude(sol, moon::MoonBody2D, [variable]; lat_idx=nothing, hours=280)

Time-Longitude Hovmöller diagram at a specific latitude (default: equator).
- `variable`: Temperature (default), Moisture, or Precipitation
"""
function plot_hovmoeller_longitude(sol, moon::MoonBody2D, ::Type{Temperature}=Temperature; lat_idx=nothing, hours=280)
    if lat_idx === nothing
        lat_idx = argmin(abs.(moon.latitudes))
    end

    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    n_times = length(sol.t) - start_idx + 1
    data_matrix = zeros(moon.n_lon, n_times)

    for (col, i) in enumerate(start_idx:length(sol.t))
        T_2d = _extract_temperature(sol, moon, i)
        data_matrix[:, col] = T_2d[lat_idx, :] .- 273.15
    end

    times_hr = sol.t[start_idx:end] ./ 3600
    lat_deg = round(moon.latitudes[lat_idx], digits=1)

    p = heatmap(times_hr, moon.longitudes, data_matrix,
        xlabel="Time (hours)",
        ylabel="Longitude (°)",
        title="Time-Longitude Temperature at Latitude $(lat_deg)°",
        color=_heatmap_color(Temperature),
        colorbar_title=_colorbar_title(Temperature),
        size=(1000, 500))

    add_eclipse_shading!(times_hr[1], times_hr[end])
    return p
end

function plot_hovmoeller_longitude(sol, moon::MoonBody2D, ::Type{Moisture}; lat_idx=nothing, hours=280)
    if lat_idx === nothing
        lat_idx = argmin(abs.(moon.latitudes))
    end

    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    n_times = length(sol.t) - start_idx + 1
    data_matrix = zeros(moon.n_lon, n_times)

    for (col, i) in enumerate(start_idx:length(sol.t))
        M_2d = _extract_moisture(sol, moon, i)
        data_matrix[:, col] = M_2d[lat_idx, :]
    end

    times_hr = sol.t[start_idx:end] ./ 3600
    lat_deg = round(moon.latitudes[lat_idx], digits=1)

    p = heatmap(times_hr, moon.longitudes, data_matrix,
        xlabel="Time (hours)",
        ylabel="Longitude (°)",
        title="Time-Longitude Moisture at Latitude $(lat_deg)°",
        color=_heatmap_color(Moisture),
        colorbar_title=_colorbar_title(Moisture),
        size=(1000, 500))

    add_eclipse_shading!(times_hr[1], times_hr[end])
    return p
end

function plot_hovmoeller_longitude(sol, moon::MoonBody2D, ::Type{Precipitation}; lat_idx=nothing, hours=280)
    if lat_idx === nothing
        lat_idx = argmin(abs.(moon.latitudes))
    end

    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    n_times = length(sol.t) - start_idx + 1
    data_matrix = zeros(moon.n_lon, n_times)

    for (col, i) in enumerate(start_idx:length(sol.t))
        T_2d = _extract_temperature(sol, moon, i)
        M_2d = _extract_moisture(sol, moon, i)
        precip = _compute_precipitation(T_2d, M_2d, moon)
        data_matrix[:, col] = precip[lat_idx, :] .* 3600
    end

    times_hr = sol.t[start_idx:end] ./ 3600
    lat_deg = round(moon.latitudes[lat_idx], digits=1)

    p = heatmap(times_hr, moon.longitudes, data_matrix,
        xlabel="Time (hours)",
        ylabel="Longitude (°)",
        title="Time-Longitude Precipitation at Latitude $(lat_deg)°",
        color=_heatmap_color(Precipitation),
        colorbar_title=_colorbar_title(Precipitation),
        size=(1000, 500))

    add_eclipse_shading!(times_hr[1], times_hr[end])
    return p
end

"""
    plot_hovmoeller_latitude(sol, moon::MoonBody2D, [variable]; lon_idx=nothing, hours=280)

Time-Latitude Hovmöller diagram at a specific longitude (default: subsolar point).
- `variable`: Temperature (default), Moisture, or Precipitation
"""
function plot_hovmoeller_latitude(sol, moon::MoonBody2D, ::Type{Temperature}=Temperature; lon_idx=nothing, hours=280)
    if lon_idx === nothing
        lon_idx = argmin(abs.(moon.longitudes))
    end

    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    n_times = length(sol.t) - start_idx + 1
    data_matrix = zeros(moon.n_lat, n_times)

    for (col, i) in enumerate(start_idx:length(sol.t))
        T_2d = _extract_temperature(sol, moon, i)
        data_matrix[:, col] = T_2d[:, lon_idx] .- 273.15
    end

    times_hr = sol.t[start_idx:end] ./ 3600
    lon_deg = round(moon.longitudes[lon_idx], digits=1)

    p = heatmap(times_hr, moon.latitudes, data_matrix,
        xlabel="Time (hours)",
        ylabel="Latitude (°)",
        title="Time-Latitude Temperature at Longitude $(lon_deg)°",
        color=_heatmap_color(Temperature),
        colorbar_title=_colorbar_title(Temperature),
        size=(1000, 500))

    add_eclipse_shading!(times_hr[1], times_hr[end])
    return p
end

function plot_hovmoeller_latitude(sol, moon::MoonBody2D, ::Type{Moisture}; lon_idx=nothing, hours=280)
    if lon_idx === nothing
        lon_idx = argmin(abs.(moon.longitudes))
    end

    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    n_times = length(sol.t) - start_idx + 1
    data_matrix = zeros(moon.n_lat, n_times)

    for (col, i) in enumerate(start_idx:length(sol.t))
        M_2d = _extract_moisture(sol, moon, i)
        data_matrix[:, col] = M_2d[:, lon_idx]
    end

    times_hr = sol.t[start_idx:end] ./ 3600
    lon_deg = round(moon.longitudes[lon_idx], digits=1)

    p = heatmap(times_hr, moon.latitudes, data_matrix,
        xlabel="Time (hours)",
        ylabel="Latitude (°)",
        title="Time-Latitude Moisture at Longitude $(lon_deg)°",
        color=_heatmap_color(Moisture),
        colorbar_title=_colorbar_title(Moisture),
        size=(1000, 500))

    add_eclipse_shading!(times_hr[1], times_hr[end])
    return p
end

function plot_hovmoeller_latitude(sol, moon::MoonBody2D, ::Type{Precipitation}; lon_idx=nothing, hours=280)
    if lon_idx === nothing
        lon_idx = argmin(abs.(moon.longitudes))
    end

    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    n_times = length(sol.t) - start_idx + 1
    data_matrix = zeros(moon.n_lat, n_times)

    for (col, i) in enumerate(start_idx:length(sol.t))
        T_2d = _extract_temperature(sol, moon, i)
        M_2d = _extract_moisture(sol, moon, i)
        precip = _compute_precipitation(T_2d, M_2d, moon)
        data_matrix[:, col] = precip[:, lon_idx] .* 3600
    end

    times_hr = sol.t[start_idx:end] ./ 3600
    lon_deg = round(moon.longitudes[lon_idx], digits=1)

    p = heatmap(times_hr, moon.latitudes, data_matrix,
        xlabel="Time (hours)",
        ylabel="Latitude (°)",
        title="Time-Latitude Precipitation at Longitude $(lon_deg)°",
        color=_heatmap_color(Precipitation),
        colorbar_title=_colorbar_title(Precipitation),
        size=(1000, 500))

    add_eclipse_shading!(times_hr[1], times_hr[end])
    return p
end

# -----------------------------------------------------------------------------
# Latitude Mean + Range Plots
# -----------------------------------------------------------------------------

"""
    plot_latitude_mean_range(sol, moon::MoonBody2D, [variable])

Plot variable by latitude with mean and min/max range.
- `variable`: Temperature (default) or Moisture
"""
function plot_latitude_mean_range(sol, moon::MoonBody2D, ::Type{Temperature}=Temperature)
    n_start = max(1, length(sol.t) - div(length(sol.t), 5))
    n_lat = moon.n_lat

    vals_mean = zeros(n_lat)
    vals_min = zeros(n_lat)
    vals_max = zeros(n_lat)

    for i in 1:n_lat
        all_vals = Float64[]
        for j in n_start:length(sol.t)
            T_2d = _extract_temperature(sol, moon, j)
            append!(all_vals, T_2d[i, :])
        end
        vals_mean[i] = mean(all_vals) - 273.15
        vals_min[i] = minimum(all_vals) - 273.15
        vals_max[i] = maximum(all_vals) - 273.15
    end

    p = plot(moon.latitudes, vals_mean,
        xlabel="Latitude (°)",
        ylabel="Temperature (°C)",
        title="Temperature by Latitude (Mean + Range)",
        linewidth=3,
        color=_plot_color(Temperature),
        label="Mean",
        legend=:topright,
        size=(800, 500))

    plot!(moon.latitudes, vals_min,
        fillrange=vals_max,
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

function plot_latitude_mean_range(sol, moon::MoonBody2D, ::Type{Moisture})
    n_start = max(1, length(sol.t) - div(length(sol.t), 5))
    n_lat = moon.n_lat

    vals_mean = zeros(n_lat)
    vals_min = zeros(n_lat)
    vals_max = zeros(n_lat)

    for i in 1:n_lat
        all_vals = Float64[]
        for j in n_start:length(sol.t)
            M_2d = _extract_moisture(sol, moon, j)
            append!(all_vals, M_2d[i, :])
        end
        vals_mean[i] = mean(all_vals)
        vals_min[i] = minimum(all_vals)
        vals_max[i] = maximum(all_vals)
    end

    p = plot(moon.latitudes, vals_mean,
        xlabel="Latitude (°)",
        ylabel="Moisture (kg/m²)",
        title="Moisture by Latitude (Mean + Range)",
        linewidth=3,
        color=_plot_color(Moisture),
        label="Mean",
        legend=:topright,
        size=(800, 500))

    plot!(moon.latitudes, vals_min,
        fillrange=vals_max,
        fillalpha=0.2,
        fillcolor=:lightblue,
        linewidth=1,
        linestyle=:dash,
        color=:lightblue,
        label="Min/Max Range")

    return p
end

# -----------------------------------------------------------------------------
# Terrain Map (no variable dispatch needed)
# -----------------------------------------------------------------------------

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

    contour!(moon.longitudes, moon.latitudes, moon.elevation,
        levels=[0.0], color=:black, linewidth=1, label="")

    return p
end