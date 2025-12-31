"""
Initial state estimation for climate simulations.
"""

# =============================================================================
# Temperature Initialization
# =============================================================================

"""
    estimate_latitude_equilibrium_temperature(lat_deg) -> Float64

Estimate equilibrium temperature based on latitude alone.
Uses radiative balance: T ∝ cos(lat)^0.25
"""
@inline function estimate_latitude_equilibrium_temperature(lat_deg::Real)
    cos_lat = max(0.1, cos(deg2rad(lat_deg)))
    return INIT_T_POLE + (INIT_T_EQUATOR - INIT_T_POLE) * cos_lat^0.25
end

"""
    estimate_altitude_adjusted_temperature(T_base, elevation) -> Float64

Adjust temperature for altitude using lapse rate.
Only positive elevations (land) contribute to cooling.
"""
@inline function estimate_altitude_adjusted_temperature(T_base::Real, elevation::Real)
    elev_meters = max(0.0, elevation) * ELEVATION_SCALE
    T_adjusted = T_base - LAPSE_RATE * elev_meters
    return max(200.0, T_adjusted)  # Floor to avoid numerical issues
end

"""
    estimate_cell_equilibrium_temperature(lat, elevation) -> Float64

Estimate equilibrium temperature for a grid cell.
Combines latitude-based temperature with altitude cooling.
"""
function estimate_cell_equilibrium_temperature(lat::Real, elevation::Real)
    T_base = estimate_latitude_equilibrium_temperature(lat)
    return estimate_altitude_adjusted_temperature(T_base, elevation)
end

"""
    initialize_equilibrium_temperature_field!(T0, moon)

Fill temperature matrix with equilibrium estimates based on latitude and elevation.
"""
function initialize_equilibrium_temperature_field!(T0::Matrix, moon)
    for i in 1:moon.n_lat
        for j in 1:moon.n_lon
            T0[i, j] = estimate_cell_equilibrium_temperature(
                moon.latitudes[i],
                moon.elevation[i, j]
            )
        end
    end
end

# =============================================================================
# Moisture Initialization
# =============================================================================

"""
    initialize_ocean_cells_to_saturation!(M0, T0, moon)

Set ocean cells to saturated moisture, land cells to zero.
"""
function initialize_ocean_cells_to_saturation!(M0::Matrix, T0::Matrix, moon)
    fill!(M0, 0.0)
    for i in 1:moon.n_lat
        for j in 1:moon.n_lon
            if moon.elevation[i, j] < 0.0
                M0[i, j] = compute_saturation_moisture_at_temperature(T0[i, j])
            end
        end
    end
end

"""
    diffuse_moisture_iteration!(M_new, M_old, moon)

Perform one iteration of moisture diffusion from ocean to land.
"""
function diffuse_moisture_iteration!(M_new::Matrix, M_old::Matrix, moon)
    n_lat, n_lon = moon.n_lat, moon.n_lon

    for i in 1:n_lat
        for j in 1:n_lon
            # Only update land cells
            if moon.elevation[i, j] >= 0.0
                total_weight = 0.0
                total_moisture = 0.0

                # North neighbor
                if i > 1
                    w = moon.moisture_transport_coeffs[i, j, 1]
                    total_weight += w
                    total_moisture += w * M_old[i-1, j]
                end

                # South neighbor
                if i < n_lat
                    w = moon.moisture_transport_coeffs[i, j, 2]
                    total_weight += w
                    total_moisture += w * M_old[i+1, j]
                end

                # East neighbor (wraps)
                j_east = mod1(j + 1, n_lon)
                w = moon.moisture_transport_coeffs[i, j, 3]
                total_weight += w
                total_moisture += w * M_old[i, j_east]

                # West neighbor (wraps)
                j_west = mod1(j - 1, n_lon)
                w = moon.moisture_transport_coeffs[i, j, 4]
                total_weight += w
                total_moisture += w * M_old[i, j_west]

                if total_weight > 0.0
                    M_new[i, j] = INIT_MOISTURE_DECAY * total_moisture / total_weight
                end
            else
                # Keep ocean cells unchanged
                M_new[i, j] = M_old[i, j]
            end
        end
    end
end

"""
    diffuse_moisture_from_ocean_to_land!(M0, moon; iterations)

Diffuse moisture from ocean cells onto land using transport coefficients.
Mountain barriers limit moisture penetration inland.
"""
function diffuse_moisture_from_ocean_to_land!(M0::Matrix, moon;
                                               iterations::Int=INIT_MOISTURE_DIFFUSE_ITERS)
    M_temp = similar(M0)

    for _ in 1:iterations
        diffuse_moisture_iteration!(M_temp, M0, moon)
        M0 .= M_temp
    end
end

"""
    cap_moisture_at_saturation_fraction!(M0, T0; fraction)

Cap moisture at a fraction of local saturation to avoid oversaturation.
"""
function cap_moisture_at_saturation_fraction!(M0::Matrix, T0::Matrix;
                                               fraction::Float64=INIT_MOISTURE_SATURATION_CAP)
    for i in eachindex(M0, T0)
        M_sat = compute_saturation_moisture_at_temperature(T0[i])
        M0[i] = min(M0[i], fraction * M_sat)
    end
end

# =============================================================================
# Combined Initialization
# =============================================================================

"""
    estimate_initial_temperature_and_moisture(moon) -> (T0::Matrix, M0::Matrix)

Estimate initial temperature and moisture fields for a MoonBody2D.

Pipeline:
1. Initialize temperature from latitude/elevation equilibrium
2. Set ocean cells to saturation moisture
3. Diffuse moisture from ocean to land
4. Cap moisture at fraction of saturation
"""
function estimate_initial_temperature_and_moisture(moon)
    n_lat, n_lon = moon.n_lat, moon.n_lon

    T0 = zeros(n_lat, n_lon)
    M0 = zeros(n_lat, n_lon)

    initialize_equilibrium_temperature_field!(T0, moon)
    initialize_ocean_cells_to_saturation!(M0, T0, moon)
    diffuse_moisture_from_ocean_to_land!(M0, moon)
    cap_moisture_at_saturation_fraction!(M0, T0)

    return T0, M0
end
