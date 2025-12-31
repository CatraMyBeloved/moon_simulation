"""
Two-layer atmosphere physics: vertical exchange and upper layer dynamics.
"""

# =============================================================================
# Ascent Conditions
# =============================================================================

"""
    is_hot_temperature_anomaly(T_surface, T_zonal_mean) -> Bool

Check if surface temperature exceeds zonal mean (hot anomaly).
"""
@inline is_hot_temperature_anomaly(T_surface::Real, T_zonal_mean::Real) = T_surface > T_zonal_mean

"""
    has_sufficient_moisture_for_deep_convection(M_surface) -> Bool

Check if there's enough surface moisture for deep convection.
"""
@inline has_sufficient_moisture_for_deep_convection(M_surface::Real) = M_surface >= MIN_MOISTURE_CONVECTION

"""
    compute_latitude_convection_factor(lat_deg) -> Float64

Calculate latitude-dependent convection factor.
Convection is easier in the tropics (cos(lat) is larger).
"""
@inline function compute_latitude_convection_factor(lat_deg::Real)
    return cos(deg2rad(lat_deg))
end

"""
    compute_upper_layer_backpressure_factor(U_above) -> Float64

Calculate suppression of ascent due to mass already aloft.
High U above suppresses further ascent.
"""
@inline function compute_upper_layer_backpressure_factor(U_above::Real)
    U_mean = 1.0
    return U_mean / (U_mean + max(U_above, U_FLOOR))
end

# =============================================================================
# Ascent Rate
# =============================================================================

"""
    compute_convective_ascent_rate(T_surf, T_zonal, U_above, M_surf, lat_deg) -> Float64

Compute vertical ascent rate from surface to upper atmosphere.

Ascent requires:
1. Surface temperature exceeding zonal mean (hot anomaly)
2. Sufficient moisture for deep convection
3. Upper atmosphere not already saturated with mass

Returns ascent rate in 1/s.
"""
@inline function compute_convective_ascent_rate(T_surf::Real, T_zonal::Real, U_above::Real,
                                                 M_surf::Real, lat_deg::Real)
    ΔT = T_surf - T_zonal

    if ΔT <= 0.0 || !has_sufficient_moisture_for_deep_convection(M_surf)
        return 0.0
    end

    lat_factor = compute_latitude_convection_factor(lat_deg)
    backpressure = compute_upper_layer_backpressure_factor(U_above)

    rate = lat_factor * (ΔT / THERMAL_RESPONSE_SCALE) * backpressure * ASCENT_RATE_MAX

    return min(rate, ASCENT_RATE_MAX)
end

# =============================================================================
# Descent Rate
# =============================================================================

"""
    compute_mass_excess_descent_rate(U_above) -> Float64

Calculate descent rate driven by excess mass aloft (U > 1).
"""
@inline function compute_mass_excess_descent_rate(U_above::Real)
    U_excess = max(0.0, U_above - 1.0)
    return MASS_DESCENT_COEFF * U_excess * U_above
end

"""
    compute_polar_sink_descent_enhancement(U_above, lat_deg) -> Float64

Calculate enhanced descent at high latitudes (polar sink).
"""
@inline function compute_polar_sink_descent_enhancement(U_above::Real, lat_deg::Real)
    if abs(lat_deg) <= POLAR_SINK_LAT
        return 0.0
    end
    polar_factor = (abs(lat_deg) - POLAR_SINK_LAT) / (90.0 - POLAR_SINK_LAT)
    return POLAR_SINK_STRENGTH * polar_factor * U_above
end

"""
    compute_total_descent_rate(U_above, lat_deg) -> Float64

Compute total vertical descent rate from upper atmosphere.

Descent is driven by:
1. Base subsidence (always present)
2. Excess mass aloft (U > 1 must descend)
3. Polar sink (enhanced descent at high latitudes)

Returns descent rate in 1/s, scaled by available mass.
"""
@inline function compute_total_descent_rate(U_above::Real, lat_deg::Real)
    rate = BASE_DESCENT_RATE
    rate += compute_mass_excess_descent_rate(U_above)
    rate += compute_polar_sink_descent_enhancement(U_above, lat_deg)
    return rate * U_above  # Scale by available mass
end

# =============================================================================
# Moisture Transfer During Vertical Motion
# =============================================================================

"""
    compute_temperature_at_lifting_condensation_level(T_surface) -> Float64

Calculate temperature where rising air reaches saturation.
"""
@inline function compute_temperature_at_lifting_condensation_level(T_surface::Real)
    return T_surface - LIFT_TEMPERATURE_DROP
end

"""
    compute_moisture_lifted_during_ascent(ascent_rate, M_surface) -> Float64

Calculate amount of moisture attempting to rise with ascending air.
"""
@inline function compute_moisture_lifted_during_ascent(ascent_rate::Real, M_surface::Real)
    return ascent_rate * M_surface * LIFT_FRACTION
end

"""
    compute_ascent_moisture_transfer(ascent_rate, M_surf, T_surf) -> (precip::Float64, M_to_upper::Float64)

Compute moisture transfer during ascent.

Rising air carries moisture upward, but most precipitates during lifting
due to adiabatic cooling. Only a fraction survives to reach the upper layer.

Returns tuple: (precipitation_during_ascent, moisture_reaching_upper_layer)
"""
@inline function compute_ascent_moisture_transfer(ascent_rate::Real, M_surf::Real, T_surf::Real)
    if ascent_rate <= 0.0 || M_surf <= 0.0
        return (0.0, 0.0)
    end

    M_lifted = compute_moisture_lifted_during_ascent(ascent_rate, M_surf)

    T_lcl = compute_temperature_at_lifting_condensation_level(T_surf)
    M_sat_lcl = compute_saturation_moisture_at_temperature(T_lcl)

    M_can_survive = M_sat_lcl * MOISTURE_SURVIVE_FRACTION
    M_to_upper = min(M_lifted, M_can_survive)
    precip_ascent = max(0.0, M_lifted - M_to_upper)

    return (precip_ascent, M_to_upper)
end

"""
    compute_descent_moisture_transfer(descent_rate, M_up, U_above) -> Float64

Compute moisture flux to surface during descent.
Descending air brings its specific moisture content back down.
"""
@inline function compute_descent_moisture_transfer(descent_rate::Real, M_up::Real, U_above::Real)
    if descent_rate <= 0.0 || M_up <= 0.0
        return 0.0
    end

    specific_M = M_up / max(U_above, U_FLOOR)
    return descent_rate * specific_M
end

# =============================================================================
# Descent Drying Effect
# =============================================================================

"""
    compute_descent_saturation_multiplier(descent_rate) -> Float64

Calculate factor by which descent suppresses surface precipitation.
Descending air warms adiabatically, raising effective saturation threshold.
Returns multiplier ≥ 1.0 for saturation moisture.
"""
@inline function compute_descent_saturation_multiplier(descent_rate::Real)
    return 1.0 + DESCENT_DRYING_SCALE * descent_rate / BASE_DESCENT_RATE
end

# =============================================================================
# Upper Mass Floor Restoration
# =============================================================================

"""
    compute_upper_mass_floor_restoration(U_current) -> Float64

Calculate restoration force when U drops below minimum floor.
Prevents numerical instability from negative or very small U.
"""
@inline function compute_upper_mass_floor_restoration(U_current::Real)
    if U_current < U_FLOOR
        return U_FLOOR_RESTORATION_RATE * (U_FLOOR - U_current)
    else
        return 0.0
    end
end
