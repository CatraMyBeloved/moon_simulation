"""
Two-Layer Atmosphere State Handling (ARCHIVED)

This file contains state vector packing/unpacking functions for the two-layer
atmosphere implementation that has been archived.

To restore: include this file after loading the main HotMoon module.
"""

# =============================================================================
# State Layout Types for Two-Layer
# =============================================================================

struct TwoLayerAtmosphereLayout <: StateLayout end
struct FullTwoLayerAtmosphereLayout <: StateLayout end

# =============================================================================
# State Vector Length for Two-Layer
# =============================================================================

"""
    state_vector_length(moon, ::TwoLayerAtmosphereLayout) -> Int

Return expected state vector length for two-layer atmosphere simulations.
"""
state_vector_length(moon, ::TwoLayerAtmosphereLayout) = 4 * moon.n_lat * moon.n_lon

"""
    state_vector_length(moon, ::FullTwoLayerAtmosphereLayout) -> Int

Return expected state vector length for full two-layer atmosphere simulations (with T_up).
"""
state_vector_length(moon, ::FullTwoLayerAtmosphereLayout) = 5 * moon.n_lat * moon.n_lon

# =============================================================================
# Upper Layer Field Unpacking
# =============================================================================

"""
    unpack_upper_mass_field(u, moon) -> Matrix view

Extract upper atmosphere mass field from two-layer state vector.
"""
@inline function unpack_upper_mass_field(u, moon)
    n_cells = moon.n_lat * moon.n_lon
    return reshape(@view(u[2n_cells+1:3n_cells]), moon.n_lat, moon.n_lon)
end

"""
    unpack_upper_moisture_field(u, moon) -> Matrix view

Extract upper atmosphere moisture field from two-layer state vector.
"""
@inline function unpack_upper_moisture_field(u, moon)
    n_cells = moon.n_lat * moon.n_lon
    return reshape(@view(u[3n_cells+1:4n_cells]), moon.n_lat, moon.n_lon)
end

"""
    unpack_upper_temperature_field(u, moon) -> Matrix view

Extract upper atmosphere temperature field from full two-layer state vector.
"""
@inline function unpack_upper_temperature_field(u, moon)
    n_cells = moon.n_lat * moon.n_lon
    return reshape(@view(u[4n_cells+1:5n_cells]), moon.n_lat, moon.n_lon)
end

# =============================================================================
# Two-Layer State Unpacking
# =============================================================================

"""
    unpack_twolayer_state(u, moon) -> NamedTuple{(:T, :M, :U, :M_up)}

Unpack two-layer atmosphere state vector into named fields.
"""
@inline function unpack_twolayer_state(u, moon)
    return (
        T = unpack_temperature_field(u, moon),
        M = unpack_moisture_field(u, moon),
        U = unpack_upper_mass_field(u, moon),
        M_up = unpack_upper_moisture_field(u, moon)
    )
end

"""
    unpack_full_twolayer_state(u, moon) -> NamedTuple{(:T, :M, :U, :M_up, :T_up)}

Unpack full two-layer atmosphere state vector into named fields (includes T_up).
"""
@inline function unpack_full_twolayer_state(u, moon)
    return (
        T = unpack_temperature_field(u, moon),
        M = unpack_moisture_field(u, moon),
        U = unpack_upper_mass_field(u, moon),
        M_up = unpack_upper_moisture_field(u, moon),
        T_up = unpack_upper_temperature_field(u, moon)
    )
end

# =============================================================================
# Two-Layer Derivative Unpacking
# =============================================================================

"""
    unpack_twolayer_derivatives(du, moon) -> NamedTuple{(:dT, :dM, :dU, :dM_up)}

Unpack derivative buffer for two-layer atmosphere simulations.
"""
@inline function unpack_twolayer_derivatives(du, moon)
    return (
        dT = unpack_temperature_field(du, moon),
        dM = unpack_moisture_field(du, moon),
        dU = unpack_upper_mass_field(du, moon),
        dM_up = unpack_upper_moisture_field(du, moon)
    )
end

"""
    unpack_full_twolayer_derivatives(du, moon) -> NamedTuple{(:dT, :dM, :dU, :dM_up, :dT_up)}

Unpack derivative buffer for full two-layer atmosphere simulations (includes dT_up).
"""
@inline function unpack_full_twolayer_derivatives(du, moon)
    return (
        dT = unpack_temperature_field(du, moon),
        dM = unpack_moisture_field(du, moon),
        dU = unpack_upper_mass_field(du, moon),
        dM_up = unpack_upper_moisture_field(du, moon),
        dT_up = unpack_upper_temperature_field(du, moon)
    )
end

# =============================================================================
# Two-Layer State Packing
# =============================================================================

"""
    pack_twolayer_state(T, M, U, M_up) -> Vector

Pack all two-layer atmosphere fields into flat state vector.
Layout: [T..., M..., U..., M_up...]
"""
function pack_twolayer_state(T::AbstractMatrix, M::AbstractMatrix,
                              U::AbstractMatrix, M_up::AbstractMatrix)
    return vcat(vec(T), vec(M), vec(U), vec(M_up))
end

"""
    pack_full_twolayer_state(T, M, U, M_up, T_up) -> Vector

Pack all full two-layer atmosphere fields into flat state vector.
Layout: [T..., M..., U..., M_up..., T_up...]
"""
function pack_full_twolayer_state(T::AbstractMatrix, M::AbstractMatrix,
                                   U::AbstractMatrix, M_up::AbstractMatrix,
                                   T_up::AbstractMatrix)
    return vcat(vec(T), vec(M), vec(U), vec(M_up), vec(T_up))
end

# =============================================================================
# Two-Layer Layout Detection
# =============================================================================

"""
    is_twolayer_state(u, moon) -> Bool

Check if state vector is from a two-layer atmosphere simulation.
"""
is_twolayer_state(u, moon) = length(u) == 4 * moon.n_lat * moon.n_lon

"""
    is_full_twolayer_state(u, moon) -> Bool

Check if state vector is from a full two-layer atmosphere simulation (with T_up).
"""
is_full_twolayer_state(u, moon) = length(u) == 5 * moon.n_lat * moon.n_lon
