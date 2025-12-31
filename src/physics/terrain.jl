"""
Terrain generation and elevation-based transport coefficients.
"""

using CoherentNoise

# =============================================================================
# Elevation Generation
# =============================================================================

"""
    create_fractal_noise_sampler(; seed, octaves, persistence) -> Sampler

Create a fractal Brownian motion noise sampler for terrain generation.
"""
function create_fractal_noise_sampler(; seed::Int=42, octaves::Int=3, persistence::Float64=0.5)
    return fbm_fractal_3d(seed=seed, octaves=octaves, persistence=persistence, lacunarity=2.0)
end

"""
    sample_terrain_elevation(lat, lon, sampler; sea_level, scale) -> Float64

Sample elevation at a geographic coordinate using fractal noise.
Uses cylindrical projection for seamless longitude wrapping.
Returns elevation where negative = ocean, positive = land.
"""
function sample_terrain_elevation(lat::Float64, lon::Float64, sampler;
                                   sea_level::Float64=0.1, scale::Float64=0.02)
    lon_rad = deg2rad(lon)

    x = lat * scale
    y = cos(lon_rad) * 90.0 * scale
    z = sin(lon_rad) * 90.0 * scale

    noise_val = sample(sampler, x, y, z)
    return noise_val - sea_level
end

# =============================================================================
# Heat Transport Coefficients
# =============================================================================

"""
    compute_slope_barrier_factor(elevation_difference) -> Float64

Calculate transport reduction from steep slopes.
Steeper slopes block transport regardless of direction.
"""
@inline function compute_slope_barrier_factor(Δelev::Real)
    return 1.0 / (1.0 + SLOPE_BARRIER_STRENGTH * Δelev^2)
end

"""
    compute_downslope_flow_bonus(elevation_difference) -> Float64

Calculate asymmetry factor for directional transport.
Positive when flowing downhill (neighbor higher), negative when uphill.
"""
@inline function compute_downslope_flow_bonus(Δelev::Real)
    return 1.0 - DOWNSLOPE_STRENGTH * tanh(Δelev * SLOPE_SCALE)
end

"""
    compute_directional_heat_transport_modifier(elev_self, elev_neighbor) -> Float64

Calculate transport coefficient modifier for heat flow FROM neighbor TO self.
Combines barrier effect and downhill preference.
"""
function compute_directional_heat_transport_modifier(elev_self::Float64, elev_neighbor::Float64)
    Δelev = elev_neighbor - elev_self
    barrier = compute_slope_barrier_factor(Δelev)
    asymmetry = compute_downslope_flow_bonus(Δelev)
    return barrier * asymmetry
end

"""
    is_ocean_cell(elevation) -> Bool

Check if a cell is ocean (elevation below sea level).
"""
@inline is_ocean_cell(elevation::Real) = elevation < 0.0

"""
    compute_ocean_transport_multiplier(elevation) -> Float64

Return transport bonus for ocean cells (currents transport heat efficiently).
"""
@inline function compute_ocean_transport_multiplier(elevation::Real)
    return is_ocean_cell(elevation) ? OCEAN_TRANSPORT_BONUS : 1.0
end

"""
    populate_heat_transport_coefficients!(coeffs, elevation, n_lat, n_lon)

Fill the 3D transport coefficient array based on elevation differences.
Directions: 1=North (i-1), 2=South (i+1), 3=East, 4=West
"""
function populate_heat_transport_coefficients!(transport_coeffs::Array{Float64,3},
                                                elevation::Matrix{Float64},
                                                n_lat::Int, n_lon::Int)
    @inbounds for i in 1:n_lat
        for j in 1:n_lon
            elev_self = elevation[i, j]
            ocean_mult = compute_ocean_transport_multiplier(elev_self)

            # North neighbor (direction 1)
            if i > 1
                elev_neighbor = elevation[i-1, j]
                transport_coeffs[i, j, 1] = ocean_mult * compute_directional_heat_transport_modifier(elev_self, elev_neighbor)
            else
                transport_coeffs[i, j, 1] = 0.0
            end

            # South neighbor (direction 2)
            if i < n_lat
                elev_neighbor = elevation[i+1, j]
                transport_coeffs[i, j, 2] = ocean_mult * compute_directional_heat_transport_modifier(elev_self, elev_neighbor)
            else
                transport_coeffs[i, j, 2] = 0.0
            end

            # East neighbor (direction 3, wraps)
            j_east = mod1(j + 1, n_lon)
            elev_neighbor = elevation[i, j_east]
            transport_coeffs[i, j, 3] = ocean_mult * compute_directional_heat_transport_modifier(elev_self, elev_neighbor)

            # West neighbor (direction 4, wraps)
            j_west = mod1(j - 1, n_lon)
            elev_neighbor = elevation[i, j_west]
            transport_coeffs[i, j, 4] = ocean_mult * compute_directional_heat_transport_modifier(elev_self, elev_neighbor)
        end
    end
end

# =============================================================================
# Moisture Transport Coefficients
# =============================================================================

"""
    compute_directional_moisture_transport_modifier(elev_self, elev_neighbor) -> Float64

Calculate transport modifier for moisture flow between cells.
Mountains block moisture more strongly than heat.
"""
function compute_directional_moisture_transport_modifier(elev_self::Float64, elev_neighbor::Float64)
    Δelev = elev_neighbor - elev_self

    barrier = 1.0 / (1.0 + MOISTURE_BARRIER_STRENGTH * Δelev^2)
    asymmetry = 1.0 - MOISTURE_DOWNSLOPE_PREFERENCE * tanh(Δelev * MOISTURE_SLOPE_SENSITIVITY)

    return barrier * asymmetry
end

"""
    populate_moisture_transport_coefficients!(coeffs, elevation, n_lat, n_lon)

Fill the 3D moisture transport coefficient array based on elevation differences.
Directions: 1=North, 2=South, 3=East, 4=West
"""
function populate_moisture_transport_coefficients!(coeffs::Array{Float64,3},
                                                    elevation::Matrix{Float64},
                                                    n_lat::Int, n_lon::Int)
    @inbounds for i in 1:n_lat
        for j in 1:n_lon
            elev_self = elevation[i, j]

            # North (direction 1)
            if i > 1
                elev_neighbor = elevation[i-1, j]
                coeffs[i, j, 1] = compute_directional_moisture_transport_modifier(elev_self, elev_neighbor)
            else
                coeffs[i, j, 1] = 0.0
            end

            # South (direction 2)
            if i < n_lat
                elev_neighbor = elevation[i+1, j]
                coeffs[i, j, 2] = compute_directional_moisture_transport_modifier(elev_self, elev_neighbor)
            else
                coeffs[i, j, 2] = 0.0
            end

            # East (direction 3, wraps)
            j_east = mod1(j + 1, n_lon)
            elev_neighbor = elevation[i, j_east]
            coeffs[i, j, 3] = compute_directional_moisture_transport_modifier(elev_self, elev_neighbor)

            # West (direction 4, wraps)
            j_west = mod1(j - 1, n_lon)
            elev_neighbor = elevation[i, j_west]
            coeffs[i, j, 4] = compute_directional_moisture_transport_modifier(elev_self, elev_neighbor)
        end
    end
end

# =============================================================================
# Temperature-Dependent Transport
# =============================================================================

"""
    compute_convection_activation(T_avg) -> Float64

Calculate convection activation factor based on average temperature.
Returns 0 when cold, 1 when warm, with smooth tanh transition.
"""
@inline @fastmath function compute_convection_activation(T_avg::Real)
    return 0.5 * (1 + tanh((T_avg - CONVECTION_THRESHOLD) / CONVECTION_WIDTH))
end

"""
    compute_heat_transport_coefficient(T_avg) -> Float64

Calculate total heat transport coefficient including convection.
Combines base conduction with temperature-activated convection.
"""
@inline @fastmath function compute_heat_transport_coefficient(T_avg::Real)
    activation = compute_convection_activation(T_avg)
    return HEAT_TRANSPORT_BASE + HEAT_TRANSPORT_CONVECTION * activation
end
