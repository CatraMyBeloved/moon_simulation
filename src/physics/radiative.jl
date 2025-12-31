"""
Radiative physics: albedo, greenhouse effect, IR transmissivity.
"""

# =============================================================================
# Ice-Albedo Feedback
# =============================================================================

"""
    compute_ice_fraction_from_temperature(T) -> Float64

Compute the fraction of surface covered by ice based on temperature.
Uses tanh transition centered at freezing point.
"""
@inline @fastmath function compute_ice_fraction_from_temperature(T::Real)
    x = (T - FREEZING_POINT) / ICE_TRANSITION_WIDTH
    return 0.5 * (1 - tanh(x))
end

"""
    compute_surface_albedo(T) -> Float64

Calculate surface albedo with ice-albedo feedback.
Smoothly interpolates between base albedo (warm) and ice albedo (cold).
"""
@inline @fastmath function compute_surface_albedo(T::Real)
    ice_frac = compute_ice_fraction_from_temperature(T)
    return BASE_ALBEDO * (1 - ice_frac) + ICE_ALBEDO * ice_frac
end

# =============================================================================
# Greenhouse Effect
# =============================================================================

"""
    compute_water_vapor_optical_depth(moisture_kg_m2) -> Float64

Calculate optical depth contribution from water vapor.
Uses logarithmic scaling (diminishing returns at high moisture).
"""
@inline function compute_water_vapor_optical_depth(moisture_kg_m2::Real)
    return WV_OPTICAL_SCALE * log(1 + max(0.0, moisture_kg_m2) / WV_REF_MOISTURE)
end

"""
    compute_total_ir_optical_depth(moisture_kg_m2) -> Float64

Calculate total infrared optical depth from permanent gases plus water vapor.
Higher τ = more absorption = stronger greenhouse effect.
"""
@inline function compute_total_ir_optical_depth(moisture_kg_m2::Real)
    τ_base = BASE_IR_OPTICAL_DEPTH
    τ_wv = compute_water_vapor_optical_depth(moisture_kg_m2)
    return τ_base + τ_wv
end

"""
    compute_total_ir_optical_depth() -> Float64

Calculate infrared optical depth for dry atmosphere (no moisture tracking).
Returns base optical depth from permanent gases only.
"""
@inline function compute_total_ir_optical_depth()
    return BASE_IR_OPTICAL_DEPTH
end

"""
    compute_ir_transmissivity_from_optical_depth(τ) -> Float64

Calculate fraction of surface IR radiation escaping to space.
Uses Eddington approximation: transmissivity = 2 / (2 + τ)
"""
@inline function compute_ir_transmissivity_from_optical_depth(τ::Real)
    return 2.0 / (2.0 + τ)
end

# =============================================================================
# Combined Radiation Calculations
# =============================================================================

"""
    compute_effective_optical_depth_at_elevation(τ_base, elevation) -> Float64

Adjust optical depth for elevation (thinner atmosphere at altitude).
"""
@inline function compute_effective_optical_depth_at_elevation(τ_base::Real, elevation::Real)
    return τ_base * (1.0 - ELEVATION_GREENHOUSE_REDUCTION * max(0.0, elevation))
end

"""
    compute_outgoing_longwave_radiation(T, transmissivity) -> Float64

Calculate outgoing longwave radiation from surface.
Q_out = ε × σ × T⁴ × transmissivity
"""
@inline @fastmath function compute_outgoing_longwave_radiation(T::Real, transmissivity::Real)
    return EMISSIVITY * STEFAN_BOLTZMANN * T^4 * transmissivity
end
