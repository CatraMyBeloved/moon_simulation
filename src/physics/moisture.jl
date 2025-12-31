"""
Moisture physics: saturation, evaporation, and precipitation.
"""

# =============================================================================
# Saturation
# =============================================================================

"""
    compute_saturation_moisture_at_temperature(T) -> Float64

Compute saturation moisture capacity using Clausius-Clapeyron relation.
Warmer air holds more moisture, colder air holds less.
At T=273K returns MOISTURE_REF_SATURATION (5.0 kg/m²).
"""
@inline @fastmath function compute_saturation_moisture_at_temperature(T::Real)
    exponent = CLAUSIUS_CLAPEYRON_SCALE * (T - MOISTURE_REF_TEMP) / MOISTURE_REF_TEMP
    return MOISTURE_REF_SATURATION * exp(clamp(exponent, -10.0, 10.0))
end

# =============================================================================
# Evaporation
# =============================================================================

"""
    can_evaporate(T, elevation) -> Bool

Check if conditions allow evaporation (ocean cell above threshold temperature).
"""
@inline function can_evaporate(T::Real, elevation::Real)
    return elevation < 0.0 && T > EVAP_THRESHOLD
end

"""
    compute_ocean_evaporation_rate(T, elevation) -> Float64

Calculate evaporation rate from ocean surface.
Only ocean cells (elevation < 0) evaporate, and only when warm enough.
Returns evaporation in kg/m²/s.
"""
@inline function compute_ocean_evaporation_rate(T::Real, elevation::Real)
    if !can_evaporate(T, elevation)
        return 0.0
    end
    return EVAP_RATE * (T - EVAP_THRESHOLD)
end

# =============================================================================
# Precipitation
# =============================================================================

"""
    compute_temperature_at_altitude(T_surface, elevation) -> Float64

Calculate temperature at altitude using lapse rate.
Only positive elevations (land) contribute to cooling.
"""
@inline function compute_temperature_at_altitude(T_surface::Real, elevation::Real)
    elev_meters = max(0.0, elevation) * ELEVATION_SCALE
    return T_surface - LAPSE_RATE * elev_meters
end

"""
    compute_orographic_precipitation_rate(moisture, T, elevation) -> Float64

Calculate precipitation rate based on orographic cooling.
When moisture exceeds saturation capacity at altitude, precipitation occurs.
Returns precipitation in kg/m²/s.
"""
@inline function compute_orographic_precipitation_rate(moisture::Real, T::Real, elevation::Real)
    if moisture <= 0.0
        return 0.0
    end

    T_altitude = compute_temperature_at_altitude(T, elevation)
    M_sat = compute_saturation_moisture_at_temperature(T_altitude)

    if moisture > M_sat
        return PRECIP_RATE * (moisture - M_sat)
    end

    return 0.0
end

# =============================================================================
# Latent Heat
# =============================================================================

"""
    compute_latent_heat_flux(precipitation, evaporation) -> Float64

Calculate net latent heat flux from phase changes.
Positive when precipitation exceeds evaporation (releases heat).
"""
@inline function compute_latent_heat_flux(precipitation::Real, evaporation::Real)
    return LATENT_HEAT * (precipitation - evaporation)
end
