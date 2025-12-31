"""
1D visualization functions for latitude-band simulations.
"""

# =============================================================================
# 1D Global Mean
# =============================================================================

function plot_1d_global_mean_temperature(sol, moon::MoonBody1D)
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

    add_day_night_shading_to_plot!(t_start_hr, t_end_hr)
    add_eclipse_shading_to_plot!(t_start_hr, t_end_hr)

    return p
end

# =============================================================================
# 1D Latitude Profile
# =============================================================================

function plot_1d_latitude_temperature_profile(sol, moon::MoonBody1D)
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

    hline!([0], linestyle=:dash, color=:blue, label="Freezing", alpha=0.5)
end

# =============================================================================
# 1D Heatmap
# =============================================================================

function plot_1d_temperature_heatmap(sol, moon::MoonBody1D)
    times_hr = sol.t ./ 3600

    temps_matrix = hcat(sol.u...)
    temps_matrix_C = temps_matrix .- 273.15

    heatmap(times_hr, moon.latitudes, temps_matrix_C,
        xlabel="Time (hours)",
        ylabel="Latitude (°)",
        title="Temperature Evolution",
        color=:thermal,
        size=(900, 400))
end

# =============================================================================
# 1D Latitude Mean + Range
# =============================================================================

function plot_1d_latitude_temperature_range(sol, moon::MoonBody1D)
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

    plot!(moon.latitudes, temps_min_C,
        fillrange=temps_max_C,
        fillalpha=0.3,
        fillcolor=:orange,
        linewidth=1,
        linestyle=:dash,
        color=:orange,
        label="Day/Night Range")

    hline!([0], linestyle=:dash, color=:blue, alpha=0.5, label="Freezing (0°C)")
    hline!([100], linestyle=:dot, color=:red, alpha=0.3, label="Boiling (100°C)")

    hspan!([0, 40], color=:green, alpha=0.1, label="")

    return p
end

# =============================================================================
# 1D Latitude Timeseries
# =============================================================================

function plot_1d_latitude_temperature_timeseries(sol, moon::MoonBody1D)
    n_points = length(sol.t)
    start_idx = max(1, n_points - 1100)

    times_hr = sol.t[start_idx:end] ./ 3600

    colors = cgrad(:thermal, moon.n_lat, categorical=true)

    p = plot(
        xlabel="Time (hours)",
        ylabel="Temperature (°C)",
        title="Temperature by Latitude Over Time",
        legend=:outerright,
        size=(1000, 500))

    for i in 1:moon.n_lat
        temps_C = [sol.u[j][i] - 273.15 for j in start_idx:n_points]
        lat_label = "$(round(Int, moon.latitudes[i]))°"
        plot!(times_hr, temps_C,
            color=colors[i],
            linewidth=1.5,
            label=lat_label,
            alpha=0.8)
    end

    t_start_hr = times_hr[1]
    t_end_hr = times_hr[end]
    add_eclipse_shading_to_plot!(t_start_hr, t_end_hr)

    hline!([0], linestyle=:dash, color=:blue, alpha=0.5, label="")

    return p
end

# =============================================================================
# 1D Summary (save all plots)
# =============================================================================

function plot_1d_summary_to_files(sol, moon::MoonBody1D, filename="summary.png")
    base = replace(filename, r"\.(png|pdf|svg)$" => "")

    p1 = plot_1d_global_mean_temperature(sol, moon)
    savefig(p1, "$(base)_global_mean.png")
    println("Saved: $(base)_global_mean.png")

    p2 = plot_1d_latitude_temperature_profile(sol, moon)
    savefig(p2, "$(base)_latitude_profile.png")
    println("Saved: $(base)_latitude_profile.png")

    p3 = plot_1d_temperature_heatmap(sol, moon)
    savefig(p3, "$(base)_heatmap.png")
    println("Saved: $(base)_heatmap.png")

    p4 = plot_1d_latitude_temperature_range(sol, moon)
    savefig(p4, "$(base)_latitude_range.png")
    println("Saved: $(base)_latitude_range.png")

    p5 = plot_1d_latitude_temperature_timeseries(sol, moon)
    savefig(p5, "$(base)_latitude_timeseries.png")
    println("Saved: $(base)_latitude_timeseries.png")
end
