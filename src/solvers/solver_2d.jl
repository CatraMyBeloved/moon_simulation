"""
2D ODE solvers for latitude-longitude climate models.

Provides two solver variants:
- Temperature-only: derivatives_temperature_only!
- Temperature + Moisture: derivatives_temperature_moisture!

Note: Two-layer atmosphere solvers have been archived to src/archived/twolayer/
"""

using DifferentialEquations

# =============================================================================
# Transport Computation Helpers
# =============================================================================

"""
    accumulate_heat_transport_from_neighbor!(transport, T, moon, i, j, ni, nj, direction)

Add heat transport contribution from a single neighbor to cell (i,j).
"""
@inline function accumulate_heat_transport_from_neighbor!(transport, T, moon, i, j, ni, nj, direction)
    T_self = T[i, j]
    T_neighbor = T[ni, nj]
    k = compute_heat_transport_coefficient(0.5 * (T_self + T_neighbor))
    k_dir = k * moon.transport_coeffs[i, j, direction]
    transport[i, j] += k_dir * (T_neighbor - T_self)
end

"""
    accumulate_heat_and_moisture_transport_from_neighbor!(heat_transport, moist_transport, T, M, moon, i, j, ni, nj, direction)

Add heat and moisture transport contributions from a single neighbor.
"""
@inline function accumulate_heat_and_moisture_transport_from_neighbor!(
        heat_transport, moist_transport, T, M, moon, i, j, ni, nj, direction)

    T_self, M_self = T[i, j], M[i, j]
    T_nb, M_nb = T[ni, nj], M[ni, nj]

    # Heat transport
    k = compute_heat_transport_coefficient(0.5 * (T_self + T_nb))
    k_dir = k * moon.transport_coeffs[i, j, direction]
    heat_transport[i, j] += k_dir * (T_nb - T_self)

    # Moisture transport
    k_moist_dir = MOISTURE_DIFFUSION * moon.moisture_transport_coeffs[i, j, direction]
    moist_transport[i, j] += k_moist_dir * (M_nb - M_self)
end

"""
    compute_cell_heat_transport!(transport, T, moon, i, j)

Accumulate heat transport from all neighbors for cell (i,j).
"""
@inline function compute_cell_heat_transport!(transport, T, moon, i, j)
    n_lat, n_lon = moon.n_lat, moon.n_lon

    # North
    if i > 1
        accumulate_heat_transport_from_neighbor!(transport, T, moon, i, j, i-1, j, DIRECTION_NORTH)
    end
    # South
    if i < n_lat
        accumulate_heat_transport_from_neighbor!(transport, T, moon, i, j, i+1, j, DIRECTION_SOUTH)
    end
    # East (wraps)
    accumulate_heat_transport_from_neighbor!(transport, T, moon, i, j, i, mod1(j+1, n_lon), DIRECTION_EAST)
    # West (wraps)
    accumulate_heat_transport_from_neighbor!(transport, T, moon, i, j, i, mod1(j-1, n_lon), DIRECTION_WEST)
end

"""
    compute_cell_heat_and_moisture_transport!(heat_transport, moist_transport, T, M, moon, i, j)

Accumulate heat and moisture transport from all neighbors for cell (i,j).
"""
@inline function compute_cell_heat_and_moisture_transport!(heat_transport, moist_transport, T, M, moon, i, j)
    n_lat, n_lon = moon.n_lat, moon.n_lon

    # North
    if i > 1
        accumulate_heat_and_moisture_transport_from_neighbor!(
            heat_transport, moist_transport, T, M, moon, i, j, i-1, j, DIRECTION_NORTH)
    end
    # South
    if i < n_lat
        accumulate_heat_and_moisture_transport_from_neighbor!(
            heat_transport, moist_transport, T, M, moon, i, j, i+1, j, DIRECTION_SOUTH)
    end
    # East (wraps)
    accumulate_heat_and_moisture_transport_from_neighbor!(
        heat_transport, moist_transport, T, M, moon, i, j, i, mod1(j+1, n_lon), DIRECTION_EAST)
    # West (wraps)
    accumulate_heat_and_moisture_transport_from_neighbor!(
        heat_transport, moist_transport, T, M, moon, i, j, i, mod1(j-1, n_lon), DIRECTION_WEST)
end

# =============================================================================
# Temperature-Only Solver
# =============================================================================

"""
    derivatives_temperature_only!(dT, temps, moon::MoonBody2D, t)

Calculate temperature derivatives for each grid cell (temperature-only simulation).
"""
function derivatives_temperature_only!(dT, temps, moon::MoonBody2D, t)
    n_lat, n_lon = moon.n_lat, moon.n_lon

    T = reshape(temps, n_lat, n_lon)
    dT_2d = reshape(dT, n_lat, n_lon)
    transport = reshape(moon._transport_cache, n_lat, n_lon)

    fill!(transport, 0.0)

    # Transport loop
    foreach_cell(n_lat, n_lon) do i, j
        @inbounds compute_cell_heat_transport!(transport, T, moon, i, j)
    end

    # Balance equations
    foreach_cell(n_lat, n_lon) do i, j
        @inbounds begin
            T_cell = T[i, j]
            lat = moon.latitudes[i]
            lon = moon.longitudes[j]
            elev = moon.elevation[i, j]

            Q_in = get_solar_2d(t, lat, lon, T_cell)

            τ_IR = compute_total_ir_optical_depth()
            τ_effective = compute_effective_optical_depth_at_elevation(τ_IR, elev)
            transmissivity = compute_ir_transmissivity_from_optical_depth(τ_effective)

            Q_out = compute_outgoing_longwave_radiation(T_cell, transmissivity)

            heat_cap = compute_biome_blended_heat_capacity(T_cell, DEFAULT_MOISTURE_ESTIMATE_KG_M2, elev)

            dT_2d[i, j] = (Q_in - Q_out + transport[i, j]) / heat_cap
        end
    end
end

# =============================================================================
# Temperature + Moisture Solver
# =============================================================================

"""
    derivatives_temperature_moisture!(du, u, moon::MoonBody2D, t)

Coupled temperature-moisture ODE system.
State vector layout: [T₁₁, ..., Tₙₘ, M₁₁, ..., Mₙₘ]
"""
function derivatives_temperature_moisture!(du, u, moon::MoonBody2D, t)
    n_lat, n_lon = moon.n_lat, moon.n_lon

    state = unpack_temperature_moisture_state(u, moon)
    derivs = unpack_temperature_moisture_derivatives(du, moon)

    T, M = state.T, state.M
    dT, dM = derivs.dT, derivs.dM

    heat_transport = reshape(moon._transport_cache, n_lat, n_lon)
    moist_transport = reshape(moon._moisture_cache, n_lat, n_lon)

    fill!(heat_transport, 0.0)
    fill!(moist_transport, 0.0)

    # Transport loop
    foreach_cell(n_lat, n_lon) do i, j
        @inbounds compute_cell_heat_and_moisture_transport!(heat_transport, moist_transport, T, M, moon, i, j)
    end

    # Balance equations
    foreach_cell(n_lat, n_lon) do i, j
        @inbounds begin
            T_cell = T[i, j]
            M_cell = max(0.0, M[i, j])
            lat = moon.latitudes[i]
            lon = moon.longitudes[j]
            elev = moon.elevation[i, j]

            Q_in = get_solar_2d(t, lat, lon, T_cell)

            τ_IR = compute_total_ir_optical_depth(M_cell)
            τ_effective = compute_effective_optical_depth_at_elevation(τ_IR, elev)
            transmissivity = compute_ir_transmissivity_from_optical_depth(τ_effective)

            Q_out = compute_outgoing_longwave_radiation(T_cell, transmissivity)

            heat_cap = compute_biome_blended_heat_capacity(T_cell, M_cell, elev)

            evap = compute_ocean_evaporation_rate(T_cell, elev)
            precip = compute_orographic_precipitation_rate(M_cell, T_cell, elev)

            Q_latent = compute_latent_heat_flux(precip, evap)

            dT[i, j] = (Q_in - Q_out + heat_transport[i, j] + Q_latent) / heat_cap
            dM[i, j] = evap - precip + moist_transport[i, j]
        end
    end
end

# =============================================================================
# Simulation Runners
# =============================================================================

"""
    run_simulation(moon::MoonBody2D, hours, T0; callback=nothing)

Run a 2D temperature-only climate simulation.
"""
function run_simulation(moon::MoonBody2D, hours::Real, T0::AbstractMatrix{<:Real}; callback=nothing)
    tspan = (0.0, hours * 3600)
    T0_flat = vec(T0)

    prob = ODEProblem(derivatives_temperature_only!, T0_flat, tspan, moon)

    if callback === nothing
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1800.0)
    else
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1800.0, callback=callback)
    end

    return sol
end

"""
    run_simulation_with_moisture(moon::MoonBody2D, hours, T0, M0; callback=nothing)

Run a coupled temperature-moisture simulation.
"""
function run_simulation_with_moisture(moon::MoonBody2D, hours::Real,
                                       T0::AbstractMatrix{<:Real},
                                       M0::AbstractMatrix{<:Real};
                                       callback=nothing)
    u0 = pack_temperature_moisture_state(T0, M0)
    tspan = (0.0, hours * 3600.0)

    prob = ODEProblem(derivatives_temperature_moisture!, u0, tspan, moon)

    if callback === nothing
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8, saveat=1800.0)
    else
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8, saveat=1800.0, callback=callback)
    end

    return sol
end

