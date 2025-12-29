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

# =============================================================================
# Animation Functions
# =============================================================================

"""
    _get_subsolar_point(t::Real) -> (lon, lat)

Calculate the subsolar point (where sun is directly overhead) at time t.
Returns (longitude, latitude) in degrees.
"""
function _get_subsolar_point(t::Real)
    # Subsolar longitude moves west as moon rotates
    subsolar_lon = mod(-360.0 * (t / ROTATION_PERIOD), 360.0)
    # Shift to -180 to 180 range
    if subsolar_lon > 180
        subsolar_lon -= 360
    end
    subsolar_lat = 0.0  # No axial tilt in this model
    return (subsolar_lon, subsolar_lat)
end

"""
    _compute_animation_clims(sol, moon::MoonBody2D, ::Type{V}, frame_indices) where V <: PlotVariable

Pre-compute fixed color limits across all animation frames for consistent visualization.
Returns (vmin, vmax) tuple.
"""
function _compute_animation_clims(sol, moon::MoonBody2D, ::Type{Temperature}, frame_indices)
    all_vals = Float64[]
    for idx in frame_indices
        T_2d = _extract_temperature(sol, moon, idx) .- 273.15
        append!(all_vals, vec(T_2d))
    end
    sorted = sort(all_vals)
    n = length(sorted)
    # Use 2nd to 98th percentile for temperature
    return (sorted[max(1, div(2*n, 100))], sorted[max(1, div(98*n, 100))])
end

function _compute_animation_clims(sol, moon::MoonBody2D, ::Type{Moisture}, frame_indices)
    max_val = 0.0
    for idx in frame_indices
        M_2d = _extract_moisture(sol, moon, idx)
        max_val = max(max_val, maximum(M_2d))
    end
    return (0.0, max_val * 1.05)
end

function _compute_animation_clims(sol, moon::MoonBody2D, ::Type{Precipitation}, frame_indices)
    all_vals = Float64[]
    for idx in frame_indices
        T_2d = _extract_temperature(sol, moon, idx)
        M_2d = _extract_moisture(sol, moon, idx)
        precip = _compute_precipitation(T_2d, M_2d, moon) .* 3600
        append!(all_vals, vec(precip))
    end
    sorted = sort(all_vals)
    n = length(sorted)
    # Use 95th percentile for precipitation (can have spikes)
    return (0.0, sorted[max(1, div(95*n, 100))])
end

"""
    animate_variable(sol, moon::MoonBody2D, ::Type{V}, filename::String; kwargs...) where V <: PlotVariable

Create an animated GIF of a single variable over time.

# Arguments
- `sol`: Solution from coupled simulation
- `moon`: MoonBody2D structure
- `V`: Variable type (Temperature, Moisture, or Precipitation)
- `filename`: Output GIF path

# Keyword Arguments
- `hours`: Number of hours to animate (default: 448 = 1 metacycle)
- `frame_skip`: Time steps between frames (default: 4 = 2 hours)
- `fps`: Frames per second in output (default: 10)
- `show_subsolar`: Show subsolar point marker (default: true)
"""
function animate_variable(sol, moon::MoonBody2D, ::Type{Temperature}, filename::String;
                          hours::Real=448, frame_skip::Int=4, fps::Int=10,
                          show_subsolar::Bool=true)
    # Calculate frame indices
    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    frame_indices = collect(start_idx:frame_skip:length(sol.t))
    n_frames = length(frame_indices)

    println("  Creating $(n_frames) frames for Temperature...")

    # Pre-compute color limits
    clims = _compute_animation_clims(sol, moon, Temperature, frame_indices)

    # Create animation
    anim = @animate for (frame_num, idx) in enumerate(frame_indices)
        t = sol.t[idx]
        time_hr = round(t / 3600, digits=1)
        eclipsed = is_eclipsed(t)

        T_2d = _extract_temperature(sol, moon, idx) .- 273.15

        eclipse_str = eclipsed ? " [ECLIPSE]" : ""
        title = "Temperature (t=$(time_hr)hr)$(eclipse_str)"

        p = heatmap(moon.longitudes, moon.latitudes, T_2d,
            xlabel="Longitude (°)",
            ylabel="Latitude (°)",
            title=title,
            color=:inferno,
            clims=clims,
            colorbar_title="°C",
            size=(800, 400),
            titlefontsize=10)

        # Add coastline
        contour!(moon.longitudes, moon.latitudes, moon.elevation,
            levels=[0.0], color=:white, linewidth=0.5, label="")

        # Add subsolar marker (yellow dot) unless eclipsed
        if show_subsolar && !eclipsed
            subsolar_lon, subsolar_lat = _get_subsolar_point(t)
            scatter!([subsolar_lon], [subsolar_lat],
                color=:yellow, markersize=8, markerstrokecolor=:orange,
                markerstrokewidth=2, label="")
        end

        # Progress
        if frame_num % 20 == 0
            print("    Frame $frame_num/$n_frames\r")
        end
    end

    println("    Encoding GIF...")
    gif(anim, filename, fps=fps)
    println("    Saved: $filename")

    return filename
end

function animate_variable(sol, moon::MoonBody2D, ::Type{Moisture}, filename::String;
                          hours::Real=448, frame_skip::Int=4, fps::Int=10,
                          show_subsolar::Bool=true)
    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    frame_indices = collect(start_idx:frame_skip:length(sol.t))
    n_frames = length(frame_indices)

    println("  Creating $(n_frames) frames for Moisture...")

    clims = _compute_animation_clims(sol, moon, Moisture, frame_indices)

    anim = @animate for (frame_num, idx) in enumerate(frame_indices)
        t = sol.t[idx]
        time_hr = round(t / 3600, digits=1)
        eclipsed = is_eclipsed(t)

        M_2d = _extract_moisture(sol, moon, idx)

        eclipse_str = eclipsed ? " [ECLIPSE]" : ""
        title = "Moisture (t=$(time_hr)hr)$(eclipse_str)"

        p = heatmap(moon.longitudes, moon.latitudes, M_2d,
            xlabel="Longitude (°)",
            ylabel="Latitude (°)",
            title=title,
            color=:Blues,
            clims=clims,
            colorbar_title="kg/m²",
            size=(800, 400),
            titlefontsize=10)

        contour!(moon.longitudes, moon.latitudes, moon.elevation,
            levels=[0.0], color=:black, linewidth=0.5, label="")

        if show_subsolar && !eclipsed
            subsolar_lon, subsolar_lat = _get_subsolar_point(t)
            scatter!([subsolar_lon], [subsolar_lat],
                color=:yellow, markersize=6, markerstrokecolor=:orange,
                markerstrokewidth=2, label="")
        end

        if frame_num % 20 == 0
            print("    Frame $frame_num/$n_frames\r")
        end
    end

    println("    Encoding GIF...")
    gif(anim, filename, fps=fps)
    println("    Saved: $filename")

    return filename
end

function animate_variable(sol, moon::MoonBody2D, ::Type{Precipitation}, filename::String;
                          hours::Real=448, frame_skip::Int=4, fps::Int=10,
                          show_subsolar::Bool=true)
    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    frame_indices = collect(start_idx:frame_skip:length(sol.t))
    n_frames = length(frame_indices)

    println("  Creating $(n_frames) frames for Precipitation...")

    clims = _compute_animation_clims(sol, moon, Precipitation, frame_indices)

    anim = @animate for (frame_num, idx) in enumerate(frame_indices)
        t = sol.t[idx]
        time_hr = round(t / 3600, digits=1)
        eclipsed = is_eclipsed(t)

        T_2d = _extract_temperature(sol, moon, idx)
        M_2d = _extract_moisture(sol, moon, idx)
        precip = _compute_precipitation(T_2d, M_2d, moon) .* 3600

        eclipse_str = eclipsed ? " [ECLIPSE]" : ""
        title = "Precipitation (t=$(time_hr)hr)$(eclipse_str)"

        p = heatmap(moon.longitudes, moon.latitudes, precip,
            xlabel="Longitude (°)",
            ylabel="Latitude (°)",
            title=title,
            color=:YlGnBu,
            clims=clims,
            colorbar_title="mm/hr",
            size=(800, 400),
            titlefontsize=10)

        contour!(moon.longitudes, moon.latitudes, moon.elevation,
            levels=[0.0], color=:black, linewidth=0.5, label="")

        if show_subsolar && !eclipsed
            subsolar_lon, subsolar_lat = _get_subsolar_point(t)
            scatter!([subsolar_lon], [subsolar_lat],
                color=:yellow, markersize=6, markerstrokecolor=:orange,
                markerstrokewidth=2, label="")
        end

        if frame_num % 20 == 0
            print("    Frame $frame_num/$n_frames\r")
        end
    end

    println("    Encoding GIF...")
    gif(anim, filename, fps=fps)
    println("    Saved: $filename")

    return filename
end

"""
    animate_combined(sol, moon::MoonBody2D, filename::String; kwargs...)

Create an animated GIF with Temperature, Moisture, and Precipitation side by side.

# Keyword Arguments
- `hours`: Number of hours to animate (default: 448 = 1 metacycle)
- `frame_skip`: Time steps between frames (default: 4 = 2 hours)
- `fps`: Frames per second (default: 10)
- `show_subsolar`: Show subsolar point marker (default: true)
"""
function animate_combined(sol, moon::MoonBody2D, filename::String;
                          hours::Real=448, frame_skip::Int=4, fps::Int=10,
                          show_subsolar::Bool=true)
    # Calculate frame indices
    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    start_idx = isnothing(start_idx) ? 1 : start_idx

    frame_indices = collect(start_idx:frame_skip:length(sol.t))
    n_frames = length(frame_indices)

    println("  Creating $(n_frames) combined frames...")

    # Pre-compute color limits for all variables
    clims_T = _compute_animation_clims(sol, moon, Temperature, frame_indices)
    clims_M = _compute_animation_clims(sol, moon, Moisture, frame_indices)
    clims_P = _compute_animation_clims(sol, moon, Precipitation, frame_indices)

    anim = @animate for (frame_num, idx) in enumerate(frame_indices)
        t = sol.t[idx]
        time_hr = round(t / 3600, digits=1)
        eclipsed = is_eclipsed(t)
        eclipse_str = eclipsed ? " [ECLIPSE]" : ""

        # Extract data
        T_2d = _extract_temperature(sol, moon, idx) .- 273.15
        M_2d = _extract_moisture(sol, moon, idx)
        T_K = T_2d .+ 273.15  # Need Kelvin for precipitation calculation
        precip = _compute_precipitation(T_K, M_2d, moon) .* 3600

        # Get subsolar point
        subsolar_lon, subsolar_lat = _get_subsolar_point(t)

        # Temperature panel
        p1 = heatmap(moon.longitudes, moon.latitudes, T_2d,
            title="Temperature (°C)",
            color=:inferno,
            clims=clims_T,
            colorbar_title="°C",
            titlefontsize=9,
            xlabel="", ylabel="Lat (°)")
        contour!(moon.longitudes, moon.latitudes, moon.elevation,
            levels=[0.0], color=:white, linewidth=0.5, label="")
        if show_subsolar && !eclipsed
            scatter!([subsolar_lon], [subsolar_lat],
                color=:yellow, markersize=6, markerstrokecolor=:orange,
                markerstrokewidth=1, label="")
        end

        # Moisture panel
        p2 = heatmap(moon.longitudes, moon.latitudes, M_2d,
            title="Moisture (kg/m²)",
            color=:Blues,
            clims=clims_M,
            colorbar_title="kg/m²",
            titlefontsize=9,
            xlabel="Lon (°)", ylabel="")
        contour!(moon.longitudes, moon.latitudes, moon.elevation,
            levels=[0.0], color=:black, linewidth=0.5, label="")
        if show_subsolar && !eclipsed
            scatter!([subsolar_lon], [subsolar_lat],
                color=:yellow, markersize=6, markerstrokecolor=:orange,
                markerstrokewidth=1, label="")
        end

        # Precipitation panel
        p3 = heatmap(moon.longitudes, moon.latitudes, precip,
            title="Precipitation (mm/hr)",
            color=:YlGnBu,
            clims=clims_P,
            colorbar_title="mm/hr",
            titlefontsize=9,
            xlabel="", ylabel="")
        contour!(moon.longitudes, moon.latitudes, moon.elevation,
            levels=[0.0], color=:black, linewidth=0.5, label="")
        if show_subsolar && !eclipsed
            scatter!([subsolar_lon], [subsolar_lat],
                color=:yellow, markersize=6, markerstrokecolor=:orange,
                markerstrokewidth=1, label="")
        end

        # Combine with overall title
        plot(p1, p2, p3,
            layout=(1, 3),
            size=(2100, 450),
            plot_title="t = $(time_hr) hr$(eclipse_str)",
            plot_titlefontsize=11)

        # Progress
        if frame_num % 20 == 0
            print("    Frame $frame_num/$n_frames\r")
        end
    end

    println("    Encoding GIF...")
    gif(anim, filename, fps=fps)
    println("    Saved: $filename")

    return filename
end

"""
    animate_all(sol, moon::MoonBody2D, output_dir::String; kwargs...)

Create all animation outputs: combined GIF plus individual variable GIFs.

# Keyword Arguments
- `hours`: Number of hours to animate (default: 448 = 1 metacycle)
- `frame_skip`: Time steps between frames (default: 4 = 2 hours)
- `fps`: Frames per second (default: 10)
"""
function animate_all(sol, moon::MoonBody2D, output_dir::String;
                     hours::Real=448, frame_skip::Int=4, fps::Int=10)

    mkpath(output_dir)

    println("Creating animations ($(hours) hours, skip=$(frame_skip))...")

    # Combined animation
    animate_combined(sol, moon, joinpath(output_dir, "climate_animation.gif");
                     hours=hours, frame_skip=frame_skip, fps=fps)

    # Individual animations
    animate_variable(sol, moon, Temperature,
                     joinpath(output_dir, "temperature_animation.gif");
                     hours=hours, frame_skip=frame_skip, fps=fps)

    animate_variable(sol, moon, Moisture,
                     joinpath(output_dir, "moisture_animation.gif");
                     hours=hours, frame_skip=frame_skip, fps=fps)

    animate_variable(sol, moon, Precipitation,
                     joinpath(output_dir, "precipitation_animation.gif");
                     hours=hours, frame_skip=frame_skip, fps=fps)

    println("All animations saved to $output_dir")
end

# =============================================================================
# Biome Visualization
# =============================================================================

"Biome classification variable (computed from T, M, elevation)"
struct Biome <: PlotVariable end

# Biome colors - 12 distinct colors for categorical visualization
# Using RGB from ColorSchemes or defining directly
const BIOME_COLORS_RGB = [
    (0.12, 0.47, 0.71),   # 0: Ocean - blue
    (1.0, 1.0, 1.0),      # 1: Ice Sheet - white
    (0.65, 0.81, 0.89),   # 2: Tundra - light blue
    (0.2, 0.63, 0.17),    # 3: Boreal Forest - dark green
    (0.76, 0.65, 0.52),   # 4: Cold Steppe - tan
    (0.13, 0.55, 0.13),   # 5: Temperate Forest - forest green
    (0.80, 0.86, 0.22),   # 6: Grassland - yellow-green
    (0.0, 0.39, 0.0),     # 7: Tropical Forest - dark green
    (1.0, 0.55, 0.0),     # 8: Savanna - orange
    (0.96, 0.87, 0.70),   # 9: Hot Desert - sand
    (0.59, 0.59, 0.59),   # 10: Mountain - gray
    (0.25, 0.60, 0.60),   # 11: Wetland - teal
]

"""
    _compute_biome_map(T_2d, M_2d, moon::MoonBody2D) -> Matrix{Int}

Compute biome classification for each cell from temperature, moisture, and elevation.
Returns matrix of biome IDs (0-11).
"""
function _compute_biome_map(T_2d, M_2d, moon::MoonBody2D)
    biome_map = zeros(Int, moon.n_lat, moon.n_lon)
    for i in 1:moon.n_lat
        for j in 1:moon.n_lon
            biome_map[i, j] = classify_biome(T_2d[i, j], M_2d[i, j], moon.elevation[i, j])
        end
    end
    return biome_map
end

"""
    plot_snapshot(sol, moon::MoonBody2D, ::Type{Biome}; time_idx=-1)

Create a biome classification map at a specific time.
Uses categorical colors with a legend showing biome names.
"""
function plot_snapshot(sol, moon::MoonBody2D, ::Type{Biome}; time_idx=-1)
    idx = time_idx < 0 ? length(sol.t) : time_idx

    T_2d = _extract_temperature(sol, moon, idx)
    M_2d = _extract_moisture(sol, moon, idx)
    biome_map = _compute_biome_map(T_2d, M_2d, moon)

    time_hr = round(sol.t[idx] / 3600, digits=1)

    # Create custom colormap from biome colors
    # Convert to a format suitable for heatmap
    biome_cmap = cgrad([RGB(c...) for c in BIOME_COLORS_RGB], NUM_BIOMES, categorical=true)

    # Plot with discrete color levels
    p = heatmap(moon.longitudes, moon.latitudes, biome_map,
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
        title="Biome Classification (t = $(time_hr) hr)",
        color=biome_cmap,
        clims=(-0.5, NUM_BIOMES - 0.5),
        colorbar=false,  # We'll add a legend instead
        size=(1200, 500))

    # Add coastline
    contour!(moon.longitudes, moon.latitudes, moon.elevation,
        levels=[0.0], color=:black, linewidth=1, label="")

    return p
end

"""
    print_biome_statistics(T_2d, M_2d, moon::MoonBody2D)

Print statistics about biome coverage.
Shows percentage of surface area covered by each biome type.
"""
function print_biome_statistics(T_2d, M_2d, moon::MoonBody2D)
    biome_map = _compute_biome_map(T_2d, M_2d, moon)

    # Calculate area-weighted coverage for each biome
    biome_areas = zeros(NUM_BIOMES)
    total_area = sum(moon.cell_areas)

    for i in 1:moon.n_lat
        for j in 1:moon.n_lon
            biome_id = biome_map[i, j]
            biome_areas[biome_id + 1] += moon.cell_areas[i, j]
        end
    end

    println("\nBiome Coverage (area-weighted):")
    println("=" ^ 40)

    for (id, name) in enumerate(BIOME_NAMES)
        pct = 100 * biome_areas[id] / total_area
        if pct > 0.1  # Only show biomes with >0.1% coverage
            bar_len = round(Int, pct / 2)  # Scale for display
            bar = "█" ^ bar_len
            println(@sprintf("  %18s: %5.1f%% %s", name, pct, bar))
        end
    end
    println("=" ^ 40)
end

"""
    plot_biome_with_legend(sol, moon::MoonBody2D; time_idx=-1)

Create a biome map with a separate legend panel showing all biome types.
Returns a combined plot with map and legend.
"""
function plot_biome_with_legend(sol, moon::MoonBody2D; time_idx=-1)
    idx = time_idx < 0 ? length(sol.t) : time_idx

    T_2d = _extract_temperature(sol, moon, idx)
    M_2d = _extract_moisture(sol, moon, idx)
    biome_map = _compute_biome_map(T_2d, M_2d, moon)

    time_hr = round(sol.t[idx] / 3600, digits=1)

    # Create custom colormap
    biome_cmap = cgrad([RGB(c...) for c in BIOME_COLORS_RGB], NUM_BIOMES, categorical=true)

    # Main biome map
    p1 = heatmap(moon.longitudes, moon.latitudes, biome_map,
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
        title="Biome Classification (t = $(time_hr) hr)",
        color=biome_cmap,
        clims=(-0.5, NUM_BIOMES - 0.5),
        colorbar=false,
        size=(900, 450))

    contour!(moon.longitudes, moon.latitudes, moon.elevation,
        levels=[0.0], color=:black, linewidth=1, label="")

    # Create legend as scatter plot with colored markers
    p2 = plot(legend=false, grid=false, axis=false, ticks=false,
        xlims=(0, 1), ylims=(0, NUM_BIOMES + 1),
        size=(250, 450))

    for (id, name) in enumerate(BIOME_NAMES)
        y_pos = NUM_BIOMES - id + 1
        color = RGB(BIOME_COLORS_RGB[id]...)

        # Color box
        scatter!([0.1], [y_pos], marker=:square, markersize=12,
            color=color, markerstrokecolor=:black, markerstrokewidth=1)

        # Label
        annotate!(0.2, y_pos, text(name, :left, 8))
    end

    # Combine plots
    return plot(p1, p2, layout=@layout([a{0.8w} b]), size=(1150, 450))
end