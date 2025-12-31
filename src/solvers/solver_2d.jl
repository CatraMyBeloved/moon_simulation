"""
2D ODE solvers for latitude-longitude climate models.

Provides three solver variants:
- Temperature-only: derivatives_temperature_only!
- Temperature + Moisture: derivatives_temperature_moisture!
- Two-layer atmosphere: derivatives_twolayer_atmosphere!
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
# Two-Layer Atmosphere Solver
# =============================================================================

"""
    compute_zonal_mean_temperature!(zonal_T, T, n_lat, n_lon)

Compute zonal (longitude) mean temperature for each latitude band.
"""
function compute_zonal_mean_temperature!(zonal_T::Vector, T, n_lat::Int, n_lon::Int)
    @inbounds for i in 1:n_lat
        total = 0.0
        for j in 1:n_lon
            total += T[i, j]
        end
        zonal_T[i] = total / n_lon
    end
end

"""
    compute_upper_layer_transport!(U_transport, M_up_transport, U, M_up, moon, i, j)

Compute upper layer mass and moisture transport for cell (i,j) using upwind scheme.
"""
@inline function compute_upper_layer_transport!(U_transport, M_up_transport, U, M_up, moon, i, j)
    n_lat, n_lon = moon.n_lat, moon.n_lon
    lat = moon.latitudes[i]

    U_here = max(U[i, j], U_FLOOR)
    M_up_here = max(M_up[i, j], 0.0)
    specific_M = M_up_here / U_here

    # Helper for upwind moisture transport
    apply_upwind_moisture = (flow, U_nb, M_up_nb) -> begin
        if flow > 0  # inflow
            return flow * max(M_up_nb, 0.0) / max(U_nb, U_FLOOR)
        else  # outflow
            return flow * specific_M
        end
    end

    # North neighbor
    if i > 1
        U_nb = max(U[i-1, j], U_FLOOR)
        ΔU = U_nb - U_here
        flow = UPPER_MERIDIONAL_COEFF * ΔU
        U_transport[i, j] += flow
        M_up_transport[i, j] += apply_upwind_moisture(flow, U_nb, M_up[i-1, j])
    end

    # South neighbor
    if i < n_lat
        U_nb = max(U[i+1, j], U_FLOOR)
        ΔU = U_nb - U_here
        flow = UPPER_MERIDIONAL_COEFF * ΔU
        U_transport[i, j] += flow
        M_up_transport[i, j] += apply_upwind_moisture(flow, U_nb, M_up[i+1, j])
    end

    # East neighbor (includes westerly bias)
    j_east = mod1(j + 1, n_lon)
    U_nb = max(U[i, j_east], U_FLOOR)
    ΔU = U_nb - U_here
    flow_pressure = UPPER_ZONAL_COEFF * ΔU
    flow_westerly = WESTERLY_BIAS_STRENGTH * abs(sin(deg2rad(lat))) * U_here
    flow = flow_pressure - flow_westerly  # westerly: outflow to east
    U_transport[i, j] += flow
    M_up_transport[i, j] += apply_upwind_moisture(flow, U_nb, M_up[i, j_east])

    # West neighbor (includes westerly bias)
    j_west = mod1(j - 1, n_lon)
    U_nb = max(U[i, j_west], U_FLOOR)
    ΔU = U_nb - U_here
    flow = UPPER_ZONAL_COEFF * ΔU + flow_westerly  # westerly: inflow from west
    U_transport[i, j] += flow
    M_up_transport[i, j] += apply_upwind_moisture(flow, U_nb, M_up[i, j_west])
end

"""
    derivatives_twolayer_atmosphere!(du, u, moon::MoonBody2D, t)

Coupled ODE system for two-layer atmosphere model.
State vector layout: [T_surf..., M_surf..., U..., M_up...]
"""
function derivatives_twolayer_atmosphere!(du, u, moon::MoonBody2D, t)
    n_lat, n_lon = moon.n_lat, moon.n_lon

    state = unpack_twolayer_state(u, moon)
    derivs = unpack_twolayer_derivatives(du, moon)

    T, M, U, M_up = state.T, state.M, state.U, state.M_up
    dT, dM, dU, dM_up = derivs.dT, derivs.dM, derivs.dU, derivs.dM_up

    # Caches
    heat_transport = reshape(moon._transport_cache, n_lat, n_lon)
    moist_transport = reshape(moon._moisture_cache, n_lat, n_lon)
    U_transport = reshape(moon._upper_U_cache, n_lat, n_lon)
    M_up_transport = reshape(moon._upper_M_cache, n_lat, n_lon)
    zonal_T = moon._zonal_mean_cache

    fill!(heat_transport, 0.0)
    fill!(moist_transport, 0.0)
    fill!(U_transport, 0.0)
    fill!(M_up_transport, 0.0)

    # Compute zonal means
    compute_zonal_mean_temperature!(zonal_T, T, n_lat, n_lon)

    # Surface transport (sequential to avoid race conditions in upper layer)
    foreach_cell_sequential(n_lat, n_lon) do i, j
        @inbounds begin
            compute_cell_heat_and_moisture_transport!(heat_transport, moist_transport, T, M, moon, i, j)
            compute_upper_layer_transport!(U_transport, M_up_transport, U, M_up, moon, i, j)
        end
    end

    # Balance equations
    foreach_cell(n_lat, n_lon) do i, j
        @inbounds begin
            T_cell = T[i, j]
            M_cell = max(0.0, M[i, j])
            U_raw = U[i, j]
            U_cell = max(U_FLOOR, U_raw)
            M_up_cell = max(0.0, M_up[i, j])
            lat = moon.latitudes[i]
            lon = moon.longitudes[j]
            elev = moon.elevation[i, j]

            # Vertical exchange
            ascent = compute_convective_ascent_rate(T_cell, zonal_T[i], U_cell, M_cell, lat)
            descent = compute_total_descent_rate(U_cell, lat)
            floor_restore = compute_upper_mass_floor_restoration(U_raw)

            dU[i, j] = ascent - descent + U_transport[i, j] + floor_restore

            # Moisture exchange
            precip_ascent, M_to_upper = compute_ascent_moisture_transfer(ascent, M_cell, T_cell)
            M_from_upper = compute_descent_moisture_transfer(descent, M_up_cell, U_cell)

            dM_up[i, j] = M_to_upper - M_from_upper + M_up_transport[i, j]

            # Surface radiation
            Q_in = get_solar_2d(t, lat, lon, T_cell)

            τ_IR = compute_total_ir_optical_depth(M_cell + M_up_cell)
            τ_effective = compute_effective_optical_depth_at_elevation(τ_IR, elev)
            transmissivity = compute_ir_transmissivity_from_optical_depth(τ_effective)

            Q_out = compute_outgoing_longwave_radiation(T_cell, transmissivity)

            heat_cap = compute_biome_blended_heat_capacity(T_cell, M_cell, elev)

            # Surface moisture with descent suppression
            evap = compute_ocean_evaporation_rate(T_cell, elev)

            drying_factor = compute_descent_saturation_multiplier(descent)
            M_sat_effective = compute_saturation_moisture_at_temperature(T_cell) * drying_factor

            precip_surf = M_cell > M_sat_effective ? PRECIP_RATE * (M_cell - M_sat_effective) : 0.0

            total_precip = precip_surf + precip_ascent
            Q_latent = compute_latent_heat_flux(total_precip, evap)

            M_lifted = ascent * M_cell * LIFT_FRACTION
            dM[i, j] = evap - precip_surf - M_lifted + M_from_upper + moist_transport[i, j]

            dT[i, j] = (Q_in - Q_out + heat_transport[i, j] + Q_latent) / heat_cap
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

"""
    run_simulation_with_twolayer_atmosphere(moon::MoonBody2D, hours, T0, M0; U0=nothing, M_up0=nothing, callback=nothing)

Run a two-layer atmosphere simulation.
"""
function run_simulation_with_twolayer_atmosphere(moon::MoonBody2D, hours::Real,
                                                  T0::AbstractMatrix{<:Real},
                                                  M0::AbstractMatrix{<:Real};
                                                  U0::Union{AbstractMatrix{<:Real},Nothing}=nothing,
                                                  M_up0::Union{AbstractMatrix{<:Real},Nothing}=nothing,
                                                  callback=nothing)
    n_lat, n_lon = moon.n_lat, moon.n_lon

    if U0 === nothing
        U0 = fill(U_INITIAL, n_lat, n_lon)
    end
    if M_up0 === nothing
        M_up0 = fill(M_UP_INITIAL, n_lat, n_lon)
    end

    u0 = pack_twolayer_state(T0, M0, U0, M_up0)
    tspan = (0.0, hours * 3600.0)

    prob = ODEProblem(derivatives_twolayer_atmosphere!, u0, tspan, moon)

    if callback === nothing
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1800.0)
    else
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1800.0, callback=callback)
    end

    return sol
end
