"""
Visualization helper functions: field extraction, solution detection, plot utilities.
"""

using Plots
using Statistics

# =============================================================================
# Plot Variable Types
# =============================================================================

abstract type PlotVariable end

struct Temperature <: PlotVariable end
struct Moisture <: PlotVariable end
struct Precipitation <: PlotVariable end
struct UpperMass <: PlotVariable end
struct UpperMoisture <: PlotVariable end
struct UpperTemperature <: PlotVariable end
struct Biome <: PlotVariable end

# =============================================================================
# Plot Configuration by Variable Type
# =============================================================================

plot_line_color(::Type{Temperature}) = :red
plot_line_color(::Type{Moisture}) = :blue
plot_line_color(::Type{Precipitation}) = :teal
plot_line_color(::Type{UpperMass}) = :purple
plot_line_color(::Type{UpperMoisture}) = :cyan
plot_line_color(::Type{UpperTemperature}) = :orange

heatmap_colorscheme(::Type{Temperature}) = :thermal
heatmap_colorscheme(::Type{Moisture}) = :Blues
heatmap_colorscheme(::Type{UpperMass}) = :viridis
heatmap_colorscheme(::Type{UpperMoisture}) = :YlGnBu
heatmap_colorscheme(::Type{UpperTemperature}) = :thermal
heatmap_colorscheme(::Type{Precipitation}) = :YlGnBu

axis_label(::Type{Temperature}) = "Temperature (°C)"
axis_label(::Type{Moisture}) = "Moisture (kg/m²)"
axis_label(::Type{Precipitation}) = "Precipitation (mm/hr)"
axis_label(::Type{UpperMass}) = "Upper Mass (relative)"
axis_label(::Type{UpperMoisture}) = "Upper Moisture (kg/m²)"
axis_label(::Type{UpperTemperature}) = "Upper Temperature (K)"

colorbar_label(::Type{Temperature}) = "°C"
colorbar_label(::Type{Moisture}) = "kg/m²"
colorbar_label(::Type{Precipitation}) = "mm/hr"
colorbar_label(::Type{UpperMass}) = "relative"
colorbar_label(::Type{UpperMoisture}) = "kg/m²"
colorbar_label(::Type{UpperTemperature}) = "K"

variable_name(::Type{Temperature}) = "Temperature"
variable_name(::Type{Moisture}) = "Moisture"
variable_name(::Type{Precipitation}) = "Precipitation"
variable_name(::Type{UpperMass}) = "Upper Mass"
variable_name(::Type{UpperMoisture}) = "Upper Moisture"
variable_name(::Type{UpperTemperature}) = "Upper Temperature"
variable_name(::Type{Biome}) = "Biome"

# Biome colors - 12 distinct colors for categorical visualization (RGB tuples 0-1)
const BIOME_COLORS_RGB = [
    (0.1, 0.2, 0.6),   # 0: Ocean (dark blue)
    (0.9, 0.95, 1.0),  # 1: Ice Sheet (white-blue)
    (0.7, 0.8, 0.7),   # 2: Tundra (gray-green)
    (0.2, 0.4, 0.3),   # 3: Boreal Forest (dark green)
    (0.8, 0.75, 0.5),  # 4: Cold Steppe (tan)
    (0.3, 0.6, 0.3),   # 5: Temperate Forest (green)
    (0.6, 0.8, 0.4),   # 6: Grassland (light green)
    (0.1, 0.5, 0.2),   # 7: Tropical Forest (dark green)
    (0.7, 0.6, 0.3),   # 8: Savanna (yellow-brown)
    (0.9, 0.8, 0.6),   # 9: Hot Desert (sand)
    (0.5, 0.45, 0.4),  # 10: Mountain (gray-brown)
    (0.3, 0.5, 0.6),   # 11: Wetland (blue-green)
]

# =============================================================================
# Solution Type Detection
# =============================================================================

is_coupled_temperature_moisture_solution(sol, moon::MoonBody2D) = length(sol.u[1]) >= 2 * moon.n_lat * moon.n_lon

is_twolayer_atmosphere_solution(sol, moon::MoonBody2D) = length(sol.u[1]) == 4 * moon.n_lat * moon.n_lon

is_full_twolayer_atmosphere_solution(sol, moon::MoonBody2D) = length(sol.u[1]) == 5 * moon.n_lat * moon.n_lon

# =============================================================================
# Field Extraction
# =============================================================================

function extract_temperature_field_celsius(sol, moon::MoonBody2D, idx::Int)
    n_cells = moon.n_lat * moon.n_lon
    if is_coupled_temperature_moisture_solution(sol, moon)
        return reshape(sol.u[idx][1:n_cells], moon.n_lat, moon.n_lon) .- 273.15
    else
        return reshape(sol.u[idx], moon.n_lat, moon.n_lon) .- 273.15
    end
end

function extract_moisture_field(sol, moon::MoonBody2D, idx::Int)
    n_cells = moon.n_lat * moon.n_lon
    return reshape(sol.u[idx][n_cells+1:2n_cells], moon.n_lat, moon.n_lon)
end

function extract_upper_mass_field(sol, moon::MoonBody2D, idx::Int)
    n_cells = moon.n_lat * moon.n_lon
    return reshape(sol.u[idx][2n_cells+1:3n_cells], moon.n_lat, moon.n_lon)
end

function extract_upper_moisture_field(sol, moon::MoonBody2D, idx::Int)
    n_cells = moon.n_lat * moon.n_lon
    return reshape(sol.u[idx][3n_cells+1:4n_cells], moon.n_lat, moon.n_lon)
end

function extract_upper_temperature_field(sol, moon::MoonBody2D, idx::Int)
    n_cells = moon.n_lat * moon.n_lon
    return reshape(sol.u[idx][4n_cells+1:5n_cells], moon.n_lat, moon.n_lon)
end

function compute_precipitation_field_mm_per_hour(sol, moon::MoonBody2D, idx::Int)
    T_2d = extract_temperature_field_celsius(sol, moon, idx) .+ 273.15
    M_2d = extract_moisture_field(sol, moon, idx)

    precip = zeros(moon.n_lat, moon.n_lon)
    for i in 1:moon.n_lat
        for j in 1:moon.n_lon
            precip[i, j] = compute_orographic_precipitation_rate(M_2d[i, j], T_2d[i, j], moon.elevation[i, j])
        end
    end
    return precip .* 3600  # Convert to mm/hr
end

# Generic field extraction by variable type
extract_field_for_plotting(sol, moon::MoonBody2D, ::Type{Temperature}, idx) = extract_temperature_field_celsius(sol, moon, idx)
extract_field_for_plotting(sol, moon::MoonBody2D, ::Type{Moisture}, idx) = extract_moisture_field(sol, moon, idx)
extract_field_for_plotting(sol, moon::MoonBody2D, ::Type{UpperMass}, idx) = extract_upper_mass_field(sol, moon, idx)
extract_field_for_plotting(sol, moon::MoonBody2D, ::Type{UpperMoisture}, idx) = extract_upper_moisture_field(sol, moon, idx)
extract_field_for_plotting(sol, moon::MoonBody2D, ::Type{UpperTemperature}, idx) = extract_upper_temperature_field(sol, moon, idx)
extract_field_for_plotting(sol, moon::MoonBody2D, ::Type{Precipitation}, idx) = compute_precipitation_field_mm_per_hour(sol, moon, idx)

# =============================================================================
# Time Range Helpers
# =============================================================================

function find_start_index_for_hours_before_end(sol, hours)
    t_end = sol.t[end]
    t_start = max(0.0, t_end - hours * 3600)
    start_idx = findfirst(t -> t >= t_start, sol.t)
    return isnothing(start_idx) ? 1 : start_idx
end

convert_time_to_hours(sol, start_idx) = sol.t[start_idx:end] ./ 3600

# =============================================================================
# Plot Decorations
# =============================================================================

function add_eclipse_shading_to_plot!(t_start_hr, t_end_hr)
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

function add_day_night_shading_to_plot!(t_start_hr, t_end_hr)
    rotation_period_hr = ROTATION_PERIOD / 3600

    for cycle_start in 0:rotation_period_hr:(t_end_hr + rotation_period_hr)
        night_start = cycle_start + 0.25 * rotation_period_hr
        night_end = cycle_start + 0.75 * rotation_period_hr

        if night_end >= t_start_hr && night_start <= t_end_hr
            vspan!([night_start, night_end], color=:navy, alpha=0.08, label="")
        end

        day1_start = cycle_start
        day1_end = cycle_start + 0.25 * rotation_period_hr
        if day1_end >= t_start_hr && day1_start <= t_end_hr
            vspan!([day1_start, day1_end], color=:yellow, alpha=0.05, label="")
        end

        day2_start = cycle_start + 0.75 * rotation_period_hr
        day2_end = cycle_start + rotation_period_hr
        if day2_end >= t_start_hr && day2_start <= t_end_hr
            vspan!([day2_start, day2_end], color=:yellow, alpha=0.05, label="")
        end
    end
end

function add_coastline_contour_to_plot!(moon::MoonBody2D)
    contour!(moon.longitudes, moon.latitudes, moon.elevation,
        levels=[0.0], color=:black, linewidth=1, label="")
end
