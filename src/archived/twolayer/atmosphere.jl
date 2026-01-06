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

"""
    compute_moisture_ascent_multiplier(M_surf) -> Float64

Calculate moisture boost factor for convective ascent.
More surface moisture = stronger convection (log scale).

At M_surf = 5 kg/m² (reference): multiplier = 1.0
At M_surf = 50 kg/m² (ocean): multiplier ≈ 2.15
"""
@inline function compute_moisture_ascent_multiplier(M_surf::Real)
    M_ref = 5.0  # reference moisture for multiplier = 1.0
    if M_surf <= M_ref
        return 1.0
    end
    return 1.0 + MOISTURE_ASCENT_SENSITIVITY * log(M_surf / M_ref)
end

# =============================================================================
# Ascent Rate
# =============================================================================

"""
    compute_convective_ascent_rate(T_surf, U_above, M_surf, lat_deg) -> Float64

Compute vertical ascent rate from surface to upper atmosphere.

Ascent requires:
1. Surface temperature exceeding absolute threshold (CONVECTION_TEMP_THRESHOLD)
2. Sufficient moisture for deep convection
3. Upper atmosphere not already saturated with mass

Higher moisture boosts ascent rate (moisture_factor), enabling ocean convection.

Returns ascent rate in 1/s.
"""
@inline function compute_convective_ascent_rate(T_surf::Real, U_above::Real,
                                                 M_surf::Real, lat_deg::Real)
    # Absolute temperature threshold (no zonal comparison)
    if T_surf < CONVECTION_TEMP_THRESHOLD || !has_sufficient_moisture_for_deep_convection(M_surf)
        return 0.0
    end

    ΔT = T_surf - CONVECTION_TEMP_THRESHOLD
    lat_factor = compute_latitude_convection_factor(lat_deg)
    backpressure = compute_upper_layer_backpressure_factor(U_above)
    moisture_factor = compute_moisture_ascent_multiplier(M_surf)

    rate = lat_factor * (ΔT / THERMAL_RESPONSE_SCALE) * moisture_factor * backpressure * ASCENT_RATE_MAX

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

# =============================================================================
# Upper Layer Temperature Physics (Phase 2 - Full Two-Layer Atmosphere)
# =============================================================================

"""
    compute_ascent_heat_transfer(ascent_rate, T_surf, T_up, U_above) -> Float64

Calculate heat transfer from surface to upper layer during ascent.
Ascending air carries sensible heat upward.

Returns temperature change rate for upper layer (K/s).
"""
@inline function compute_ascent_heat_transfer(ascent_rate::Real, T_surf::Real,
                                               T_up::Real, U_above::Real)
    if ascent_rate <= 0.0
        return 0.0
    end

    ΔT = T_surf - T_up
    heat_flux = ascent_rate * ΔT * ASCENT_HEAT_TRANSFER_FRACTION
    # Normalize by upper layer mass to get temperature change
    return heat_flux / max(U_above, U_FLOOR)
end

"""
    compute_upper_layer_latent_heating(precip_ascent, U_above) -> Float64

Calculate latent heat release in upper layer from condensation during ascent.
Moisture that precipitates during ascent releases latent heat to the upper layer.

Returns temperature change rate for upper layer (K/s).
"""
@inline function compute_upper_layer_latent_heating(precip_ascent::Real, U_above::Real)
    if precip_ascent <= 0.0
        return 0.0
    end

    # Latent heat released (J/m²/s)
    Q = LATENT_HEAT * precip_ascent
    # Convert to temperature change rate using upper layer heat capacity
    return Q / (max(U_above, U_FLOOR) * UPPER_LAYER_HEAT_CAPACITY)
end

"""
    compute_upper_layer_radiative_cooling(T_up) -> Float64

Calculate radiative cooling of upper layer toward equilibrium temperature.
Upper layer loses heat to space via infrared radiation.

Returns temperature change rate for upper layer (K/s), typically negative.
"""
@inline function compute_upper_layer_radiative_cooling(T_up::Real)
    ΔT = T_up - T_UP_EQUILIBRIUM
    return -UPPER_RADIATIVE_COOLING_RATE * ΔT
end

"""
    compute_descent_heat_transfer(descent_rate, T_up, U_above) -> Float64

Calculate heat loss from upper layer during descent.
Descending air carries heat back to surface (adiabatically warms).

Returns temperature change rate for upper layer (K/s), typically negative.
"""
@inline function compute_descent_heat_transfer(descent_rate::Real, T_up::Real, U_above::Real)
    if descent_rate <= 0.0
        return 0.0
    end

    # Descending air warms adiabatically - this heat comes from the upper layer
    # The descended air would be warmer than T_up by DESCENT_ADIABATIC_WARMING
    # This represents heat leaving the upper layer
    normalized_descent = descent_rate / BASE_DESCENT_RATE
    heat_loss_rate = descent_rate * DESCENT_ADIABATIC_WARMING * normalized_descent
    return -heat_loss_rate / max(U_above, U_FLOOR)
end

"""
    compute_vertical_instability_factor(T_surf, T_up) -> Float64

Calculate convective instability factor based on vertical temperature gradient.
Hot surface with cold upper layer = unstable = enhanced convection.

Returns multiplication factor for ascent rate (≥ 1.0).
"""
@inline function compute_vertical_instability_factor(T_surf::Real, T_up::Real)
    ΔT_vertical = T_surf - T_up
    instability = ΔT_vertical / VERTICAL_INSTABILITY_SCALE
    return 1.0 + max(0.0, instability)
end

"""
    compute_temperature_descent_factor(T_up, T_up_zonal_mean) -> Float64

Calculate descent enhancement due to cold upper layer anomaly.
Cold upper air is denser and sinks faster.

Returns multiplication factor for descent rate (≥ 1.0).
"""
@inline function compute_temperature_descent_factor(T_up::Real, T_up_zonal_mean::Real)
    cold_anomaly = max(0.0, T_up_zonal_mean - T_up)
    return 1.0 + T_DESCENT_SENSITIVITY * cold_anomaly / 10.0
end

"""
    clamp_upper_temperature(T_up) -> Float64

Clamp upper temperature to valid range for numerical stability.
"""
@inline function clamp_upper_temperature(T_up::Real)
    return clamp(T_up, T_UP_FLOOR, T_UP_CEILING)
end
