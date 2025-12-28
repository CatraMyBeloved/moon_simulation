"""
Physics functions: albedo, greenhouse effect, feedbacks, material properties
"""

"""
    get_heat_capacity(lat_deg::Float64)

Return heat capacity (J⋅m⁻²⋅K⁻¹) for a given latitude.

The heat capacity varies by zone to reflect different surface types:
- Equatorial (0-30°): Ocean/wetland with high thermal inertia
- Mid-latitude (30-60°): Mixed terrain with moderate thermal inertia
- Polar (60-90°): Dry ice/rock with low thermal inertia

This creates realistic thermal response differences: equatorial regions
change temperature slowly (buffered by water), while polar regions
respond quickly to changes in solar input.
"""
function get_heat_capacity(lat_deg::Real)
    abs_lat = abs(lat_deg)
    if abs_lat < ZONE_EQUATORIAL_END
        return HEAT_CAPACITY_EQUATORIAL
    elseif abs_lat < ZONE_MIDLAT_END
        return HEAT_CAPACITY_MIDLAT
    else
        return HEAT_CAPACITY_POLAR
    end
end

"""
    get_heat_capacity(lat_deg, lon_deg, elevation)

Return heat capacity based on elevation and latitude.

Terrain types by elevation:
- Ocean (elev < 0): High thermal inertia from water
- Wetlands (0 to ~0.15): Moist soil, high capacity
- Normal terrain (~0.15 to ~0.5): Latitude-dependent (equator vs polar)
- Mountains (> ~0.5): Bare rock, low capacity (fast response)

Smooth transitions between zones using tanh blending.
"""
function get_heat_capacity(lat_deg::Real, lon_deg::Real, elevation::Real)
    # Ocean
    if elevation < 0.0
        # Deeper = more thermal mass, but cap it
        depth_factor = clamp(1.0 - elevation, 1.0, 2.0)
        return HEAT_CAPACITY_OCEAN * depth_factor
    end

    # Get base land capacity from latitude
    land_cap = get_heat_capacity(lat_deg)

    # Wetlands zone (0 to ELEV_WETLAND_END) - blend from wetland to normal
    if elevation < ELEV_WETLAND_END
        t = elevation / ELEV_WETLAND_END  # 0 at coast, 1 at wetland end
        return HEAT_CAPACITY_WETLAND * (1 - t) + land_cap * t
    end

    # Mountain zone (above ELEV_MOUNTAIN_START) - blend from normal to rock
    if elevation > ELEV_MOUNTAIN_START
        # How far into mountain zone (0 at start, approaches 1 at high peaks)
        t = (elevation - ELEV_MOUNTAIN_START) / (1.0 - ELEV_MOUNTAIN_START)
        t = clamp(t, 0.0, 1.0)
        return land_cap * (1 - t) + HEAT_CAPACITY_MOUNTAIN * t
    end

    # Normal terrain - latitude-dependent
    return land_cap
end

"""
    get_transport_coefficient(T_avg)

Calculate heat transport coefficient with convection activation.

At low temperatures, only conduction occurs (slow transport).
At higher temperatures, convection kicks in (fast transport).
Uses tanh for smooth transition, consistent with other feedbacks.

# Arguments
- `T_avg`: Average temperature between two cells (Kelvin)

# Returns
- Transport coefficient (W⋅m⁻²⋅K⁻¹)
"""
function get_transport_coefficient(T_avg::Real)
    # Convection activation: 0 when cold, 1 when warm
    activation = 0.5 * (1 + tanh((T_avg - CONVECTION_THRESHOLD) / CONVECTION_WIDTH))

    # Total transport = base conduction + temperature-dependent convection
    return HEAT_TRANSPORT_BASE + HEAT_TRANSPORT_CONVECTION * activation
end

"""
    get_albedo(T)

Calculate surface albedo with ice-albedo feedback.

Uses a tanh transition centered at freezing point to smoothly interpolate
between ice albedo (cold) and base albedo (warm).

# Arguments
- `T`: Temperature in Kelvin

# Returns
- Albedo (0-1)
"""
function get_albedo(T::Real)
    x = (T - FREEZING_POINT) / ICE_TRANSITION_WIDTH
    ice_frac = 0.5 * (1 - tanh(x))
    return BASE_ALBEDO * (1 - ice_frac) + ICE_ALBEDO * ice_frac
end

"""
    get_greenhouse(T)

Calculate greenhouse effect with water vapor feedback.

Uses a tanh function to model increased water vapor at higher temperatures,
which strengthens the greenhouse effect.

# Arguments
- `T`: Temperature in Kelvin

# Returns
- Greenhouse factor (0-1), higher = more heat trapped
"""
function get_greenhouse(T::Real)
    temperature_delta = T - WATER_VAPOR_REF_TEMP
    water_vapor_boost = WATER_VAPOR_STRENGTH * tanh(temperature_delta / WATER_VAPOR_SCALE)
    return clamp(BASE_GREENHOUSE + water_vapor_boost, 0.2, 0.85)
end

# =============================================================================
# Terrain Generation
# =============================================================================

"""
    generate_elevation(lat, lon, sampler; sea_level=0.1, scale=0.02)

Generate elevation at a lat/lon point using fractal noise.
Uses cylindrical projection (cos/sin of longitude) for seamless wrapping.

Returns elevation where negative = ocean, positive = land.
- `sea_level`: shifts threshold (higher = more ocean)
- `scale`: controls feature size (smaller = larger continents)
"""
function generate_elevation(lat::Float64, lon::Float64, sampler; sea_level::Float64=0.1, scale::Float64=0.02)
    # Convert to cylindrical coordinates for seamless longitude wrapping
    lon_rad = deg2rad(lon)

    # Sample noise on a cylinder: (lat, cos(lon), sin(lon))
    # Scale controls feature size - smaller scale = larger landmasses
    x = lat * scale
    y = cos(lon_rad) * 90.0 * scale
    z = sin(lon_rad) * 90.0 * scale

    # Get fractal noise value (typically in [-1, 1] range)
    noise_val = sample(sampler, x, y, z)

    # Shift by sea level to control land/ocean ratio
    return noise_val - sea_level
end

"""
    create_elevation_sampler(; seed=42, octaves=3, persistence=0.5)

Create a fractal noise sampler for elevation generation.
- `octaves`: more = finer detail, less = smoother terrain
- `persistence`: how much each octave contributes (lower = smoother)
"""
function create_elevation_sampler(; seed::Int=42, octaves::Int=3, persistence::Float64=0.5)
    return fbm_fractal_3d(seed=seed, octaves=octaves, persistence=persistence, lacunarity=2.0)
end

# =============================================================================
# Elevation-Based Transport
# =============================================================================

"""
    calculate_directional_transport(elev_self, elev_neighbor)

Calculate transport coefficient modifier for flow FROM neighbor TO self.

Includes two effects:
1. Barrier effect: Steep slopes reduce transport in both directions
2. Asymmetry: Downhill flow is easier than uphill (katabatic winds, gravity)

Returns a multiplier (typically 0.2 to 2.0) for the base transport coefficient.
"""
function calculate_directional_transport(elev_self::Float64, elev_neighbor::Float64)
    Δelev = elev_neighbor - elev_self  # positive = neighbor is higher (uphill from us)

    # Barrier effect: steep slopes block transport regardless of direction
    # Using 1/(1 + k*Δe²) so very steep slopes approach zero transport
    barrier = 1.0 / (1.0 + SLOPE_BARRIER_STRENGTH * Δelev^2)

    # Asymmetry: downhill flow enhanced, uphill flow reduced
    # When Δelev < 0 (neighbor is lower, we're looking downhill): bonus
    # When Δelev > 0 (neighbor is higher, we're looking uphill): penalty
    # tanh gives smooth transition bounded to [-1, 1]
    asymmetry = 1.0 - DOWNSLOPE_STRENGTH * tanh(Δelev * SLOPE_SCALE)

    return barrier * asymmetry
end

"""
    calculate_transport_coefficients!(transport_coeffs, elevation, n_lat, n_lon)

Fill the transport coefficient array based on elevation differences.
Directions: 1=toward lower lat index (i-1), 2=toward higher lat index (i+1), 3=East, 4=West

Also applies ocean transport bonus for underwater cells.
"""
function calculate_transport_coefficients!(transport_coeffs::Array{Float64,3},
                                           elevation::Matrix{Float64},
                                           n_lat::Int, n_lon::Int)
    for i in 1:n_lat
        for j in 1:n_lon
            elev_self = elevation[i, j]

            # Ocean bonus: water transports heat more efficiently
            ocean_mult = elev_self < 0.0 ? OCEAN_TRANSPORT_BONUS : 1.0

            # Neighbor at i-1 (lower latitude index)
            if i > 1
                elev_neighbor = elevation[i-1, j]
                transport_coeffs[i, j, 1] = ocean_mult * calculate_directional_transport(elev_self, elev_neighbor)
            else
                transport_coeffs[i, j, 1] = 0.0  # No neighbor at boundary
            end

            # Neighbor at i+1 (higher latitude index)
            if i < n_lat
                elev_neighbor = elevation[i+1, j]
                transport_coeffs[i, j, 2] = ocean_mult * calculate_directional_transport(elev_self, elev_neighbor)
            else
                transport_coeffs[i, j, 2] = 0.0  # No neighbor at boundary
            end

            # East neighbor (wraps around)
            j_east = mod1(j + 1, n_lon)
            elev_neighbor = elevation[i, j_east]
            transport_coeffs[i, j, 3] = ocean_mult * calculate_directional_transport(elev_self, elev_neighbor)

            # West neighbor (wraps around)
            j_west = mod1(j - 1, n_lon)
            elev_neighbor = elevation[i, j_west]
            transport_coeffs[i, j, 4] = ocean_mult * calculate_directional_transport(elev_self, elev_neighbor)
        end
    end
end