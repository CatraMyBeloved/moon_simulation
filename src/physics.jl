"""
Physics functions: albedo, greenhouse effect, feedbacks, material properties
"""

"""
    get_heat_capacity(lat_deg::Float64)

Return heat capacity (J⋅m⁻²⋅K⁻¹) for a given latitude.

The heat capacity varies by zone to reflect different surface types:
- Equatorial (0-30°): Hot, dry land with low thermal inertia (fast response)
- Mid-latitude (30-60°): Mixed terrain with moderate thermal inertia
- Polar (60-90°): Ice sheets and frozen ocean with high thermal inertia (slow response)

This creates realistic thermal response differences: equatorial regions
respond quickly to diurnal heating cycles, while polar ice masses
change temperature slowly due to their large thermal mass.
"""
@inline function get_heat_capacity(lat_deg::Real)
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
@inline function get_heat_capacity(lat_deg::Real, lon_deg::Real, elevation::Real)
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
@inline @fastmath function get_transport_coefficient(T_avg::Real)
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
@inline @fastmath function get_albedo(T::Real)
    x = (T - FREEZING_POINT) / ICE_TRANSITION_WIDTH
    ice_frac = 0.5 * (1 - tanh(x))
    return BASE_ALBEDO * (1 - ice_frac) + ICE_ALBEDO * ice_frac
end

"""
    get_ir_optical_depth(M)

Calculate infrared optical depth from atmospheric composition.

The optical depth determines how opaque the atmosphere is to thermal radiation.
Higher τ = more absorption = stronger greenhouse effect.

# Arguments
- `M`: Moisture content in kg/m² (water vapor contribution)

# Returns
- IR optical depth τ (dimensionless, typically 0.5-3.0)
"""
@inline function get_ir_optical_depth(M::Real)
    # Base optical depth from permanent gases (CO₂, CH₄, etc.)
    τ_base = BASE_IR_OPTICAL_DEPTH

    # Water vapor adds optical depth (log scale - diminishing returns at high M)
    # At M=5 kg/m²: adds ~0.55, at M=20 kg/m²: adds ~1.2
    τ_wv = WV_OPTICAL_SCALE * log(1 + max(0.0, M) / WV_REF_MOISTURE)

    return τ_base + τ_wv
end

"""
    get_ir_optical_depth()

Calculate infrared optical depth for dry atmosphere (no moisture tracking).
Uses base optical depth only.

# Returns
- IR optical depth τ (dimensionless)
"""
@inline function get_ir_optical_depth()
    return BASE_IR_OPTICAL_DEPTH
end

"""
    get_ir_transmissivity(τ)

Calculate fraction of surface IR radiation that escapes to space.

Uses the Eddington approximation for a gray atmosphere:
    transmissivity = 2 / (2 + τ)

This comes from solving the two-stream radiative transfer equations.

# Arguments
- `τ`: IR optical depth

# Returns
- Transmissivity (0-1), fraction of surface emission reaching space

# Examples
- τ = 0: transmissivity = 1.0 (no atmosphere, all escapes)
- τ = 1: transmissivity = 0.67 (33% blocked)
- τ = 2: transmissivity = 0.50 (50% blocked)
- τ = 4: transmissivity = 0.33 (67% blocked)
"""
@inline function get_ir_transmissivity(τ::Real)
    return 2.0 / (2.0 + τ)
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
    @inbounds for i in 1:n_lat
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

# =============================================================================
# Moisture System
# =============================================================================

"""
    get_saturation_moisture(T)

Compute saturation moisture capacity at temperature T.
Based on Clausius-Clapeyron: saturation increases exponentially with temperature.

At T=273K: returns MOISTURE_REF_SATURATION (5.0 kg/m²)
Warmer air holds more moisture, colder air holds less.

# Arguments
- `T`: Temperature in Kelvin

# Returns
- Saturation moisture (kg/m²)
"""
@inline @fastmath function get_saturation_moisture(T::Real)
    exponent = CLAUSIUS_CLAPEYRON_SCALE * (T - MOISTURE_REF_TEMP) / MOISTURE_REF_TEMP
    return MOISTURE_REF_SATURATION * exp(clamp(exponent, -10.0, 10.0))
end

"""
    get_evaporation(T, elevation)

Compute evaporation rate from surface.
Only ocean cells (elevation < 0) evaporate, and only when warm enough.

# Arguments
- `T`: Surface temperature in Kelvin
- `elevation`: Normalized elevation (negative = ocean)

# Returns
- Evaporation rate (kg/m²/s)
"""
@inline function get_evaporation(T::Real, elevation::Real)
    # Only ocean evaporates
    if elevation >= 0.0
        return 0.0
    end

    # Only above threshold temperature
    if T <= EVAP_THRESHOLD
        return 0.0
    end

    # Linear increase with temperature above threshold
    return EVAP_RATE * (T - EVAP_THRESHOLD)
end

"""
    get_precipitation(M, T, elevation)

Compute precipitation rate based on saturation.
Uses lapse rate to compute effective temperature at altitude (orographic cooling).
When moisture exceeds saturation capacity at altitude, precipitation occurs.

# Arguments
- `M`: Current moisture content (kg/m²)
- `T`: Surface temperature in Kelvin
- `elevation`: Normalized elevation

# Returns
- Precipitation rate (kg/m²/s)
"""
@inline function get_precipitation(M::Real, T::Real, elevation::Real)
    # No precipitation if no moisture
    if M <= 0.0
        return 0.0
    end

    # Compute temperature at altitude (orographic cooling)
    # Only land elevations contribute to lapse rate cooling
    elev_meters = max(0.0, elevation) * ELEVATION_SCALE
    T_altitude = T - LAPSE_RATE * elev_meters

    # Saturation at altitude temperature
    M_sat = get_saturation_moisture(T_altitude)

    # Precipitate excess moisture
    if M > M_sat
        return PRECIP_RATE * (M - M_sat)
    end

    return 0.0
end

"""
    moisture_directional_transport(elev_self, elev_neighbor)

Compute moisture transport modifier between two cells.
Mountains block moisture more strongly than they block heat.

# Arguments
- `elev_self`: Elevation of current cell
- `elev_neighbor`: Elevation of neighbor cell

# Returns
- Transport modifier (0 to ~1)
"""
function moisture_directional_transport(elev_self::Float64, elev_neighbor::Float64)
    delta_elev = elev_neighbor - elev_self

    # Barrier effect: steep slopes strongly block moisture
    barrier = 1.0 / (1.0 + MOISTURE_BARRIER_STRENGTH * delta_elev^2)

    # Slight downslope preference (moist air flows downhill)
    asymmetry = 1.0 - 0.3 * tanh(delta_elev * 3.0)

    return barrier * asymmetry
end

"""
    calculate_moisture_transport_coefficients!(coeffs, elevation, n_lat, n_lon)

Pre-compute directional moisture transport modifiers.
Mountains block moisture more strongly than heat.
Directions: 1=North, 2=South, 3=East, 4=West
"""
function calculate_moisture_transport_coefficients!(coeffs::Array{Float64,3},
                                                     elevation::Matrix{Float64},
                                                     n_lat::Int, n_lon::Int)
    @inbounds for i in 1:n_lat
        for j in 1:n_lon
            elev_self = elevation[i, j]

            # North (direction 1)
            if i > 1
                elev_neighbor = elevation[i-1, j]
                coeffs[i, j, 1] = moisture_directional_transport(elev_self, elev_neighbor)
            else
                coeffs[i, j, 1] = 0.0
            end

            # South (direction 2)
            if i < n_lat
                elev_neighbor = elevation[i+1, j]
                coeffs[i, j, 2] = moisture_directional_transport(elev_self, elev_neighbor)
            else
                coeffs[i, j, 2] = 0.0
            end

            # East (direction 3) - wraps
            j_east = mod1(j + 1, n_lon)
            elev_neighbor = elevation[i, j_east]
            coeffs[i, j, 3] = moisture_directional_transport(elev_self, elev_neighbor)

            # West (direction 4) - wraps
            j_west = mod1(j - 1, n_lon)
            elev_neighbor = elevation[i, j_west]
            coeffs[i, j, 4] = moisture_directional_transport(elev_self, elev_neighbor)
        end
    end
end

# =============================================================================
# Biome System
# =============================================================================

"""
    classify_biome(T, M, elevation)

Classify a cell into one of 12 biome types based on temperature, moisture, and elevation.

Returns biome ID (0-11):
- 0: Ocean (elevation < 0)
- 1: Ice Sheet (T < -10°C)
- 2: Tundra (-10°C ≤ T < 0°C)
- 3: Boreal Forest (0°C ≤ T < 10°C, wet)
- 4: Cold Steppe (0°C ≤ T < 10°C, dry)
- 5: Temperate Forest (10°C ≤ T < 20°C, wet)
- 6: Grassland (10°C ≤ T < 20°C, dry)
- 7: Tropical Forest (T ≥ 20°C, wet)
- 8: Savanna (T ≥ 20°C, moderate moisture)
- 9: Hot Desert (T ≥ 20°C, dry)
- 10: Mountain (elevation > 0.5)
- 11: Wetland (0 < elevation < 0.1)

# Arguments
- `T`: Temperature in Kelvin
- `M`: Moisture in kg/m²
- `elevation`: Normalized elevation (-1 to 1)

# Returns
- Biome ID (Int)
"""
function classify_biome(T::Real, M::Real, elevation::Real)
    # Elevation-based biomes take priority
    if elevation < 0.0
        return BIOME_OCEAN
    elseif elevation > BIOME_ELEV_MOUNTAIN
        return BIOME_MOUNTAIN
    elseif elevation < BIOME_ELEV_WETLAND
        return BIOME_WETLAND
    end

    # Temperature-based classification
    if T < BIOME_T_ICE
        return BIOME_ICE_SHEET
    elseif T < BIOME_T_TUNDRA
        return BIOME_TUNDRA
    elseif T < BIOME_T_COLD
        # Cold biomes: boreal forest vs cold steppe
        return M >= BIOME_M_DRY ? BIOME_BOREAL_FOREST : BIOME_COLD_STEPPE
    elseif T < BIOME_T_TEMPERATE
        # Temperate biomes: forest vs grassland
        return M >= BIOME_M_MODERATE ? BIOME_TEMPERATE_FOREST : BIOME_GRASSLAND
    else
        # Hot biomes: tropical forest vs savanna vs desert
        if M >= BIOME_M_WET
            return BIOME_TROPICAL_FOREST
        elseif M >= BIOME_M_DESERT
            return BIOME_SAVANNA
        else
            return BIOME_HOT_DESERT
        end
    end
end

"""
    get_heat_capacity_biome(T, M, elevation)

Calculate heat capacity based on climate conditions using smooth biome blending.

Instead of discrete biome boundaries, this function smoothly blends between
adjacent biome heat capacities using tanh transitions. This prevents numerical
artifacts from sharp discontinuities.

# Arguments
- `T`: Temperature in Kelvin
- `M`: Moisture in kg/m²
- `elevation`: Normalized elevation

# Returns
- Heat capacity (J⋅m⁻²⋅K⁻¹)
"""
@inline @fastmath function get_heat_capacity_biome(T::Real, M::Real, elevation::Real)
    # Ocean - handle separately with depth scaling
    if elevation < 0.0
        depth_factor = clamp(1.0 - elevation, 1.0, 2.0)
        return BIOME_HEAT_CAPACITY[BIOME_OCEAN + 1] * depth_factor
    end

    # Mountain transition (smooth blend above BIOME_ELEV_MOUNTAIN)
    mountain_weight = 0.5 * (1 + tanh((elevation - BIOME_ELEV_MOUNTAIN) / BIOME_ELEV_WIDTH))
    if mountain_weight > 0.99
        return BIOME_HEAT_CAPACITY[BIOME_MOUNTAIN + 1]
    end

    # Wetland transition (smooth blend below BIOME_ELEV_WETLAND)
    wetland_weight = 0.5 * (1 - tanh((elevation - BIOME_ELEV_WETLAND) / BIOME_ELEV_WIDTH))

    # Temperature weights for blending between climate zones
    # ice_weight: fraction in ice sheet zone
    ice_weight = 0.5 * (1 - tanh((T - BIOME_T_ICE) / BIOME_T_WIDTH))

    # tundra_weight: fraction in tundra zone (between ice and cold)
    tundra_lower = 0.5 * (1 + tanh((T - BIOME_T_ICE) / BIOME_T_WIDTH))
    tundra_upper = 0.5 * (1 - tanh((T - BIOME_T_TUNDRA) / BIOME_T_WIDTH))
    tundra_weight = tundra_lower * tundra_upper

    # cold_weight: fraction in cold zone (0-10°C)
    cold_lower = 0.5 * (1 + tanh((T - BIOME_T_TUNDRA) / BIOME_T_WIDTH))
    cold_upper = 0.5 * (1 - tanh((T - BIOME_T_COLD) / BIOME_T_WIDTH))
    cold_weight = cold_lower * cold_upper

    # temperate_weight: fraction in temperate zone (10-20°C)
    temp_lower = 0.5 * (1 + tanh((T - BIOME_T_COLD) / BIOME_T_WIDTH))
    temp_upper = 0.5 * (1 - tanh((T - BIOME_T_TEMPERATE) / BIOME_T_WIDTH))
    temperate_weight = temp_lower * temp_upper

    # hot_weight: fraction in hot zone (>20°C)
    hot_weight = 0.5 * (1 + tanh((T - BIOME_T_TEMPERATE) / BIOME_T_WIDTH))

    # Moisture weights for wet/dry classification within each zone
    # dry_weight: fraction in dry conditions
    dry_weight = 0.5 * (1 - tanh((M - BIOME_M_DRY) / BIOME_M_WIDTH))
    # moderate_weight: fraction in moderate conditions
    mod_lower = 0.5 * (1 + tanh((M - BIOME_M_DESERT) / BIOME_M_WIDTH))
    mod_upper = 0.5 * (1 - tanh((M - BIOME_M_WET) / BIOME_M_WIDTH))
    moderate_weight = mod_lower * mod_upper
    # wet_weight: fraction in wet conditions
    wet_weight = 0.5 * (1 + tanh((M - BIOME_M_MODERATE) / BIOME_M_WIDTH))

    # Calculate heat capacity for each temperature zone
    # Cold zone: blend between boreal forest and cold steppe
    hc_cold = wet_weight * BIOME_HEAT_CAPACITY[BIOME_BOREAL_FOREST + 1] +
              (1 - wet_weight) * BIOME_HEAT_CAPACITY[BIOME_COLD_STEPPE + 1]

    # Temperate zone: blend between forest and grassland
    hc_temperate = wet_weight * BIOME_HEAT_CAPACITY[BIOME_TEMPERATE_FOREST + 1] +
                   (1 - wet_weight) * BIOME_HEAT_CAPACITY[BIOME_GRASSLAND + 1]

    # Hot zone: blend between tropical forest, savanna, and desert
    desert_weight = 0.5 * (1 - tanh((M - BIOME_M_DESERT) / BIOME_M_WIDTH))
    forest_weight = 0.5 * (1 + tanh((M - BIOME_M_WET) / BIOME_M_WIDTH))
    savanna_weight = 1.0 - desert_weight - forest_weight
    savanna_weight = max(0.0, savanna_weight)  # Ensure non-negative

    hc_hot = forest_weight * BIOME_HEAT_CAPACITY[BIOME_TROPICAL_FOREST + 1] +
             savanna_weight * BIOME_HEAT_CAPACITY[BIOME_SAVANNA + 1] +
             desert_weight * BIOME_HEAT_CAPACITY[BIOME_HOT_DESERT + 1]

    # Blend across temperature zones
    hc_climate = ice_weight * BIOME_HEAT_CAPACITY[BIOME_ICE_SHEET + 1] +
                 tundra_weight * BIOME_HEAT_CAPACITY[BIOME_TUNDRA + 1] +
                 cold_weight * hc_cold +
                 temperate_weight * hc_temperate +
                 hot_weight * hc_hot

    # Apply wetland and mountain modifiers
    hc_land = wetland_weight * BIOME_HEAT_CAPACITY[BIOME_WETLAND + 1] +
              (1 - wetland_weight) * hc_climate

    # Final blend with mountain
    hc_final = mountain_weight * BIOME_HEAT_CAPACITY[BIOME_MOUNTAIN + 1] +
               (1 - mountain_weight) * hc_land

    return hc_final
end

# =============================================================================
# Initial State Estimation
# =============================================================================

"""
    estimate_equilibrium_temperature(lat, elevation)

Estimate equilibrium temperature based on latitude and elevation.

Uses simplified radiative balance: temperature scales with cos(lat)^0.25
for latitude dependence, then applies lapse rate cooling for altitude.

# Arguments
- `lat`: Latitude in degrees
- `elevation`: Normalized elevation

# Returns
- Temperature in Kelvin
"""
function estimate_equilibrium_temperature(lat::Real, elevation::Real)
    # Latitude-based temperature (radiative balance gives T ∝ cos(lat)^0.25)
    cos_lat = max(0.1, cos(deg2rad(lat)))
    T_base = INIT_T_POLE + (INIT_T_EQUATOR - INIT_T_POLE) * cos_lat^0.25

    # Altitude cooling (lapse rate)
    elev_meters = max(0.0, elevation) * ELEVATION_SCALE
    T_altitude = T_base - LAPSE_RATE * elev_meters

    return max(200.0, T_altitude)  # Floor at 200K to avoid numerical issues
end

"""
    estimate_initial_state(moon::MoonBody2D)

Estimate initial temperature and moisture fields for a MoonBody2D.

Temperature is estimated from latitude/elevation equilibrium.
Moisture is diffused from ocean cells onto land, with mountain barriers.

# Arguments
- `moon`: MoonBody2D structure with elevation data

# Returns
- Tuple (T0, M0) of initial temperature and moisture matrices
"""
function estimate_initial_state(moon)
    n_lat, n_lon = moon.n_lat, moon.n_lon

    # Step 1: Estimate equilibrium temperature
    T0 = zeros(n_lat, n_lon)
    for i in 1:n_lat
        for j in 1:n_lon
            T0[i, j] = estimate_equilibrium_temperature(
                moon.latitudes[i],
                moon.elevation[i, j]
            )
        end
    end

    # Step 2: Initialize moisture - ocean cells start saturated, land starts dry
    M0 = zeros(n_lat, n_lon)
    for i in 1:n_lat
        for j in 1:n_lon
            if moon.elevation[i, j] < 0.0
                # Ocean: start at saturation
                M0[i, j] = get_saturation_moisture(T0[i, j])
            end
        end
    end

    # Step 3: Diffuse moisture from ocean onto land
    for iter in 1:INIT_MOISTURE_DIFFUSE_ITERS
        M_new = copy(M0)
        for i in 1:n_lat
            for j in 1:n_lon
                # Only update land cells
                if moon.elevation[i, j] >= 0.0
                    # Weighted average of neighbors
                    total_weight = 0.0
                    total_moisture = 0.0

                    # North neighbor
                    if i > 1
                        w = moon.moisture_transport_coeffs[i, j, 1]
                        total_weight += w
                        total_moisture += w * M0[i-1, j]
                    end

                    # South neighbor
                    if i < n_lat
                        w = moon.moisture_transport_coeffs[i, j, 2]
                        total_weight += w
                        total_moisture += w * M0[i+1, j]
                    end

                    # East neighbor (wraps)
                    j_east = mod1(j + 1, n_lon)
                    w = moon.moisture_transport_coeffs[i, j, 3]
                    total_weight += w
                    total_moisture += w * M0[i, j_east]

                    # West neighbor (wraps)
                    j_west = mod1(j - 1, n_lon)
                    w = moon.moisture_transport_coeffs[i, j, 4]
                    total_weight += w
                    total_moisture += w * M0[i, j_west]

                    if total_weight > 0.0
                        # Decay factor creates gradient from ocean
                        M_new[i, j] = INIT_MOISTURE_DECAY * total_moisture / total_weight
                    end
                end
            end
        end
        M0 = M_new
    end

    # Step 4: Cap moisture at 80% of local saturation
    for i in 1:n_lat
        for j in 1:n_lon
            M_sat = get_saturation_moisture(T0[i, j])
            M0[i, j] = min(M0[i, j], 0.8 * M_sat)
        end
    end

    return T0, M0
end