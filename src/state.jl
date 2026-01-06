"""
State vector packing and unpacking for ODE solvers.

Provides type-safe abstractions for working with flattened state vectors
used by DifferentialEquations.jl solvers.

Note: Two-layer atmosphere state handling has been archived to src/archived/twolayer/
"""

# =============================================================================
# State Layout Types
# =============================================================================

abstract type StateLayout end

struct TemperatureOnlyLayout <: StateLayout end
struct TemperatureMoistureLayout <: StateLayout end

# =============================================================================
# State Vector Length Calculation
# =============================================================================

"""
    state_vector_length(moon, ::TemperatureOnlyLayout) -> Int

Return expected state vector length for temperature-only simulations.
"""
state_vector_length(moon, ::TemperatureOnlyLayout) = moon.n_lat * moon.n_lon

"""
    state_vector_length(moon, ::TemperatureMoistureLayout) -> Int

Return expected state vector length for coupled temperature-moisture simulations.
"""
state_vector_length(moon, ::TemperatureMoistureLayout) = 2 * moon.n_lat * moon.n_lon

# =============================================================================
# State Vector Validation
# =============================================================================

"""
    validate_state_vector_length(u, moon, layout::StateLayout)

Validate that state vector has correct length for the given layout.
Throws ArgumentError if length doesn't match.
"""
function validate_state_vector_length(u, moon, layout::StateLayout)
    expected = state_vector_length(moon, layout)
    actual = length(u)
    if actual != expected
        error("State vector length mismatch: expected $expected for $(typeof(layout)), got $actual")
    end
end

# =============================================================================
# Field Unpacking (Read Access)
# =============================================================================

"""
    unpack_temperature_field(u, moon) -> Matrix view

Extract temperature field as a reshaped view of the state vector.
Works with any layout (temperature is always first n_cells elements).
"""
@inline function unpack_temperature_field(u, moon)
    n_cells = moon.n_lat * moon.n_lon
    return reshape(@view(u[1:n_cells]), moon.n_lat, moon.n_lon)
end

"""
    unpack_moisture_field(u, moon) -> Matrix view

Extract moisture field from coupled state vector.
Assumes moisture is stored after temperature.
"""
@inline function unpack_moisture_field(u, moon)
    n_cells = moon.n_lat * moon.n_lon
    return reshape(@view(u[n_cells+1:2n_cells]), moon.n_lat, moon.n_lon)
end

# =============================================================================
# Full State Unpacking
# =============================================================================

"""
    unpack_temperature_only_state(u, moon) -> NamedTuple{(:T,)}

Unpack temperature-only state vector into named fields.
"""
@inline function unpack_temperature_only_state(u, moon)
    return (T = unpack_temperature_field(u, moon),)
end

"""
    unpack_temperature_moisture_state(u, moon) -> NamedTuple{(:T, :M)}

Unpack coupled temperature-moisture state vector into named fields.
"""
@inline function unpack_temperature_moisture_state(u, moon)
    return (
        T = unpack_temperature_field(u, moon),
        M = unpack_moisture_field(u, moon)
    )
end

# =============================================================================
# Derivative Buffer Unpacking
# =============================================================================

"""
    unpack_temperature_only_derivatives(du, moon) -> NamedTuple{(:dT,)}

Unpack derivative buffer for temperature-only simulations.
"""
@inline function unpack_temperature_only_derivatives(du, moon)
    return (dT = unpack_temperature_field(du, moon),)
end

"""
    unpack_temperature_moisture_derivatives(du, moon) -> NamedTuple{(:dT, :dM)}

Unpack derivative buffer for coupled temperature-moisture simulations.
"""
@inline function unpack_temperature_moisture_derivatives(du, moon)
    return (
        dT = unpack_temperature_field(du, moon),
        dM = unpack_moisture_field(du, moon)
    )
end

# =============================================================================
# State Vector Packing
# =============================================================================

"""
    pack_temperature_only_state(T) -> Vector

Pack temperature matrix into flat state vector.
"""
function pack_temperature_only_state(T::AbstractMatrix)
    return vec(T)
end

"""
    pack_temperature_moisture_state(T, M) -> Vector

Pack temperature and moisture matrices into flat state vector.
Layout: [T₁₁, T₂₁, ..., Tₙₘ, M₁₁, M₂₁, ..., Mₙₘ]
"""
function pack_temperature_moisture_state(T::AbstractMatrix, M::AbstractMatrix)
    return vcat(vec(T), vec(M))
end

# =============================================================================
# Layout Detection
# =============================================================================

"""
    detect_state_layout(u, moon) -> StateLayout

Detect the state layout based on state vector length.
"""
function detect_state_layout(u, moon)
    n_cells = moon.n_lat * moon.n_lon
    len = length(u)

    if len == n_cells
        return TemperatureOnlyLayout()
    elseif len == 2 * n_cells
        return TemperatureMoistureLayout()
    else
        error("Unknown state layout: length $len for $n_cells cell grid")
    end
end

"""
    is_temperature_only_state(u, moon) -> Bool

Check if state vector is from a temperature-only simulation.
"""
is_temperature_only_state(u, moon) = length(u) == moon.n_lat * moon.n_lon

"""
    is_temperature_moisture_state(u, moon) -> Bool

Check if state vector is from a coupled temperature-moisture simulation.
"""
is_temperature_moisture_state(u, moon) = length(u) == 2 * moon.n_lat * moon.n_lon
