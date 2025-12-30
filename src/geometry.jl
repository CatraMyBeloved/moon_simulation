"""
Solar geometry: zenith angles, eclipse detection
"""

"""
    get_zenith(t, lat_deg)

Calculate solar zenith angle at given time and latitude.

# Arguments
- `t`: Time in seconds
- `lat_deg`: Latitude in degrees

# Returns
- Zenith angle in radians, or `nothing` if sun is below horizon
"""
@inline function get_zenith(t::Real, lat_deg::Real)
    phase = (t / ROTATION_PERIOD) % 1.0 #Where in the day are we?
    lat = deg2rad(lat_deg) #Convert latitude to radians for calculations
    hour_angle = 2π * phase # Hour angle in radians, aka at this time of day, how high is the sun?
    cos_z = cos(lat) * cos(hour_angle) # cos(lat) provides the max height of the sun at this latitude
    if cos_z <= 0
        return nothing
    end
    return acos(cos_z) # convert to actual angle
end

"""
    is_eclipsed(t)

Check if the moon is currently in eclipse (behind the parent planet).

Eclipse occurs at orbital_phase = 0.5 (halfway through orbit).

# Arguments
- `t`: Time in seconds

# Returns
- `true` if eclipsed, `false` otherwise
"""
@inline function is_eclipsed(t::Real)
    orbital_phase = (t % ORBITAL_PERIOD) / ORBITAL_PERIOD
    eclipse_half_width = (ECLIPSE_DURATION / ORBITAL_PERIOD) / 2

    return abs(orbital_phase - 0.5) < eclipse_half_width
end

"""
    get_solar(t, lat_deg, T)

Calculate solar heating at given time, latitude, and temperature (1D version).
"""
@inline function get_solar(t::Real, lat_deg::Real, T::Real)
    if is_eclipsed(t)
        return 0.0
    end

    zenith = get_zenith(t, lat_deg)
    if zenith === nothing
        return 0.0
    end

    cos_zenith = cos(zenith)
    albedo = get_albedo(T)
    air_mass = 1.0 / (cos_zenith + 0.01)
    transmission = exp(-SOLAR_OPTICAL_DEPTH * air_mass)

    return SOLAR_CONSTANT * cos_zenith * (1 - albedo) * transmission
end

# =============================================================================
# 2D versions - account for longitude
# =============================================================================

"""
    get_zenith_2d(t, lat_deg, lon_deg)

Calculate solar zenith angle for a 2D grid cell.

Longitude adds a phase offset: each longitude experiences noon at a different time.
At t=0, lon=0° is at noon; lon=180° is at midnight.

# Returns
- Zenith angle in radians, or `nothing` if sun is below horizon
"""
@inline function get_zenith_2d(t::Real, lat_deg::Real, lon_deg::Real)
    # Hour angle = rotation phase + longitude offset
    # This is the only difference from the 1D version
    hour_angle = 2π * (t / ROTATION_PERIOD) + deg2rad(lon_deg)

    lat = deg2rad(lat_deg)
    cos_z = cos(lat) * cos(hour_angle)

    if cos_z <= 0
        return nothing  # Sun below horizon
    end
    return acos(cos_z)
end

"""
    get_solar_2d(t, lat_deg, lon_deg, T)

Calculate solar heating for a 2D grid cell.

Same physics as get_solar, but uses longitude-aware zenith calculation.
"""
@inline function get_solar_2d(t::Real, lat_deg::Real, lon_deg::Real, T::Real)
    if is_eclipsed(t)
        return 0.0
    end

    zenith = get_zenith_2d(t, lat_deg, lon_deg)
    if zenith === nothing
        return 0.0
    end

    cos_zenith = cos(zenith)
    albedo = get_albedo(T)
    air_mass = 1.0 / (cos_zenith + 0.01)
    transmission = exp(-SOLAR_OPTICAL_DEPTH * air_mass)

    return SOLAR_CONSTANT * cos_zenith * (1 - albedo) * transmission
end