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
    apply_upwind_moisture_transport(flow, U_nb, M_up_nb, specific_M_local) -> Float64

Helper for upwind moisture transport in upper layer.
Avoids closure allocation by passing specific_M as parameter.
"""
@inline function apply_upwind_moisture_transport(flow::Real, U_nb::Real, M_up_nb::Real, specific_M_local::Real)
    if flow > 0  # inflow from neighbor
        return flow * max(M_up_nb, 0.0) / max(U_nb, U_FLOOR)
    else  # outflow to neighbor
        return flow * specific_M_local
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

    # North neighbor
    if i > 1
        U_nb = max(U[i-1, j], U_FLOOR)
        ΔU = U_nb - U_here
        flow = UPPER_MERIDIONAL_COEFF * ΔU
        U_transport[i, j] += flow
        M_up_transport[i, j] += apply_upwind_moisture_transport(flow, U_nb, M_up[i-1, j], specific_M)
    end

    # South neighbor
    if i < n_lat
        U_nb = max(U[i+1, j], U_FLOOR)
        ΔU = U_nb - U_here
        flow = UPPER_MERIDIONAL_COEFF * ΔU
        U_transport[i, j] += flow
        M_up_transport[i, j] += apply_upwind_moisture_transport(flow, U_nb, M_up[i+1, j], specific_M)
    end

    # East neighbor (includes westerly bias)
    j_east = mod1(j + 1, n_lon)
    U_nb = max(U[i, j_east], U_FLOOR)
    ΔU = U_nb - U_here
    flow_pressure = UPPER_ZONAL_COEFF * ΔU
    flow_westerly = WESTERLY_BIAS_STRENGTH * abs(sin(deg2rad(lat))) * U_here
    flow = flow_pressure - flow_westerly  # westerly: outflow to east
    U_transport[i, j] += flow
    M_up_transport[i, j] += apply_upwind_moisture_transport(flow, U_nb, M_up[i, j_east], specific_M)

    # West neighbor (includes westerly bias)
    j_west = mod1(j - 1, n_lon)
    U_nb = max(U[i, j_west], U_FLOOR)
    ΔU = U_nb - U_here
    flow = UPPER_ZONAL_COEFF * ΔU + flow_westerly  # westerly: inflow from west
    U_transport[i, j] += flow
    M_up_transport[i, j] += apply_upwind_moisture_transport(flow, U_nb, M_up[i, j_west], specific_M)
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

            dM_up_raw = M_to_upper - M_from_upper + M_up_transport[i, j]
            # Prevent M_up from going negative: if already near zero and derivative is negative, clamp
            if M_up_cell <= 0.01 && dM_up_raw < 0.0
                dM_up[i, j] = max(dM_up_raw, -M_up_cell * 0.1)  # Soft floor
            else
                dM_up[i, j] = dM_up_raw
            end

            # Surface radiation
            Q_in = get_solar_2d(t, lat, lon, T_cell)

            τ_IR = compute_total_ir_optical_depth(M_cell + M_up_cell)
            τ_effective = compute_effective_optical_depth_at_elevation(τ_IR, elev)
            transmissivity = compute_ir_transmissivity_from_optical_depth(τ_effective)

            Q_out = compute_outgoing_longwave_radiation(T_cell, transmissivity)

            heat_cap = compute_biome_blended_heat_capacity(T_cell, M_cell, elev)

            # Surface moisture with descent suppression and orographic effect
            evap = compute_ocean_evaporation_rate(T_cell, elev)

            # Apply lapse rate cooling for elevated terrain (orographic effect)
            elev_meters = max(0.0, elev) * ELEVATION_SCALE
            T_at_altitude = T_cell - LAPSE_RATE * elev_meters

            drying_factor = compute_descent_saturation_multiplier(descent)
            M_sat_effective = compute_saturation_moisture_at_temperature(T_at_altitude) * drying_factor

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
# Full Two-Layer Atmosphere Solver (with T_up)
# =============================================================================

"""
    apply_upwind_temperature_transport(flow, U_nb, T_up_nb, specific_T_local) -> Float64

Helper for upwind temperature transport in upper layer.
"""
@inline function apply_upwind_temperature_transport(flow::Real, U_nb::Real, T_up_nb::Real, specific_T_local::Real)
    if flow > 0  # inflow from neighbor
        return flow * T_up_nb / max(U_nb, U_FLOOR)
    else  # outflow to neighbor
        return flow * specific_T_local
    end
end

"""
    compute_upper_layer_full_transport!(U_transport, M_up_transport, T_up_transport, U, M_up, T_up, moon, i, j)

Compute upper layer mass, moisture, and temperature transport for cell (i,j).
"""
@inline function compute_upper_layer_full_transport!(U_transport, M_up_transport, T_up_transport,
                                                      U, M_up, T_up, moon, i, j)
    n_lat, n_lon = moon.n_lat, moon.n_lon
    lat = moon.latitudes[i]

    U_here = max(U[i, j], U_FLOOR)
    M_up_here = max(M_up[i, j], 0.0)
    T_up_here = T_up[i, j]
    specific_M = M_up_here / U_here
    specific_T = T_up_here / U_here

    # North neighbor
    if i > 1
        U_nb = max(U[i-1, j], U_FLOOR)
        ΔU = U_nb - U_here
        flow = UPPER_MERIDIONAL_COEFF * ΔU
        U_transport[i, j] += flow
        M_up_transport[i, j] += apply_upwind_moisture_transport(flow, U_nb, M_up[i-1, j], specific_M)
        T_up_transport[i, j] += apply_upwind_temperature_transport(flow, U_nb, T_up[i-1, j], specific_T)
    end

    # South neighbor
    if i < n_lat
        U_nb = max(U[i+1, j], U_FLOOR)
        ΔU = U_nb - U_here
        flow = UPPER_MERIDIONAL_COEFF * ΔU
        U_transport[i, j] += flow
        M_up_transport[i, j] += apply_upwind_moisture_transport(flow, U_nb, M_up[i+1, j], specific_M)
        T_up_transport[i, j] += apply_upwind_temperature_transport(flow, U_nb, T_up[i+1, j], specific_T)
    end

    # East neighbor (includes westerly bias)
    j_east = mod1(j + 1, n_lon)
    U_nb = max(U[i, j_east], U_FLOOR)
    ΔU = U_nb - U_here
    flow_pressure = UPPER_ZONAL_COEFF * ΔU
    flow_westerly = WESTERLY_BIAS_STRENGTH * abs(sin(deg2rad(lat))) * U_here
    flow = flow_pressure - flow_westerly
    U_transport[i, j] += flow
    M_up_transport[i, j] += apply_upwind_moisture_transport(flow, U_nb, M_up[i, j_east], specific_M)
    T_up_transport[i, j] += apply_upwind_temperature_transport(flow, U_nb, T_up[i, j_east], specific_T)

    # West neighbor (includes westerly bias)
    j_west = mod1(j - 1, n_lon)
    U_nb = max(U[i, j_west], U_FLOOR)
    ΔU = U_nb - U_here
    flow = UPPER_ZONAL_COEFF * ΔU + flow_westerly
    U_transport[i, j] += flow
    M_up_transport[i, j] += apply_upwind_moisture_transport(flow, U_nb, M_up[i, j_west], specific_M)
    T_up_transport[i, j] += apply_upwind_temperature_transport(flow, U_nb, T_up[i, j_west], specific_T)
end

"""
    derivatives_full_twolayer_atmosphere!(du, u, moon::MoonBody2D, t)

Coupled ODE system for full two-layer atmosphere model (with upper temperature).
State vector layout: [T_surf..., M_surf..., U..., M_up..., T_up...]
"""
function derivatives_full_twolayer_atmosphere!(du, u, moon::MoonBody2D, t)
    n_lat, n_lon = moon.n_lat, moon.n_lon

    state = unpack_full_twolayer_state(u, moon)
    derivs = unpack_full_twolayer_derivatives(du, moon)

    T, M, U, M_up, T_up = state.T, state.M, state.U, state.M_up, state.T_up
    dT, dM, dU, dM_up, dT_up = derivs.dT, derivs.dM, derivs.dU, derivs.dM_up, derivs.dT_up

    # Caches
    heat_transport = reshape(moon._transport_cache, n_lat, n_lon)
    moist_transport = reshape(moon._moisture_cache, n_lat, n_lon)
    U_transport = reshape(moon._upper_U_cache, n_lat, n_lon)
    M_up_transport = reshape(moon._upper_M_cache, n_lat, n_lon)
    T_up_transport = reshape(moon._upper_T_cache, n_lat, n_lon)
    zonal_T = moon._zonal_mean_cache
    zonal_T_up = moon._zonal_T_up_cache

    fill!(heat_transport, 0.0)
    fill!(moist_transport, 0.0)
    fill!(U_transport, 0.0)
    fill!(M_up_transport, 0.0)
    fill!(T_up_transport, 0.0)

    # Compute zonal means for both surface and upper temperatures
    compute_zonal_mean_temperature!(zonal_T, T, n_lat, n_lon)
    compute_zonal_mean_temperature!(zonal_T_up, T_up, n_lat, n_lon)

    # Surface and upper layer transport (sequential to avoid race conditions)
    foreach_cell_sequential(n_lat, n_lon) do i, j
        @inbounds begin
            compute_cell_heat_and_moisture_transport!(heat_transport, moist_transport, T, M, moon, i, j)
            compute_upper_layer_full_transport!(U_transport, M_up_transport, T_up_transport,
                                                 U, M_up, T_up, moon, i, j)
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
            T_up_cell = clamp_upper_temperature(T_up[i, j])
            lat = moon.latitudes[i]
            lon = moon.longitudes[j]
            elev = moon.elevation[i, j]

            # Vertical exchange with temperature feedback
            # Instability factor: hot surface + cold upper = enhanced convection
            instability = compute_vertical_instability_factor(T_cell, T_up_cell)
            base_ascent = compute_convective_ascent_rate(T_cell, zonal_T[i], U_cell, M_cell, lat)
            ascent = base_ascent * instability

            # Temperature-driven descent: cold upper air sinks faster
            T_descent_factor = compute_temperature_descent_factor(T_up_cell, zonal_T_up[i])
            base_descent = compute_total_descent_rate(U_cell, lat)
            descent = base_descent * T_descent_factor

            floor_restore = compute_upper_mass_floor_restoration(U_raw)

            dU[i, j] = ascent - descent + U_transport[i, j] + floor_restore

            # Moisture exchange
            precip_ascent, M_to_upper = compute_ascent_moisture_transfer(ascent, M_cell, T_cell)
            M_from_upper = compute_descent_moisture_transfer(descent, M_up_cell, U_cell)

            dM_up_raw = M_to_upper - M_from_upper + M_up_transport[i, j]
            if M_up_cell <= 0.01 && dM_up_raw < 0.0
                dM_up[i, j] = max(dM_up_raw, -M_up_cell * 0.1)
            else
                dM_up[i, j] = dM_up_raw
            end

            # Upper temperature evolution
            # 1. Latent heating from condensation during ascent
            latent_heating = compute_upper_layer_latent_heating(precip_ascent, U_cell)

            # 2. Radiative cooling to space
            radiative_cooling = compute_upper_layer_radiative_cooling(T_up_cell)

            # 3. Heat transfer from ascent (sensible heat)
            ascent_heating = compute_ascent_heat_transfer(ascent, T_cell, T_up_cell, U_cell)

            # 4. Heat loss during descent
            descent_cooling = compute_descent_heat_transfer(descent, T_up_cell, U_cell)

            dT_up_raw = latent_heating + radiative_cooling + ascent_heating + descent_cooling + T_up_transport[i, j]
            # Clamp to prevent runaway
            if T_up_cell <= T_UP_FLOOR + 5.0 && dT_up_raw < 0.0
                dT_up[i, j] = max(dT_up_raw, -(T_up_cell - T_UP_FLOOR) * 0.1)
            elseif T_up_cell >= T_UP_CEILING - 5.0 && dT_up_raw > 0.0
                dT_up[i, j] = min(dT_up_raw, (T_UP_CEILING - T_up_cell) * 0.1)
            else
                dT_up[i, j] = dT_up_raw
            end

            # Surface radiation
            Q_in = get_solar_2d(t, lat, lon, T_cell)

            τ_IR = compute_total_ir_optical_depth(M_cell + M_up_cell)
            τ_effective = compute_effective_optical_depth_at_elevation(τ_IR, elev)
            transmissivity = compute_ir_transmissivity_from_optical_depth(τ_effective)

            Q_out = compute_outgoing_longwave_radiation(T_cell, transmissivity)

            heat_cap = compute_biome_blended_heat_capacity(T_cell, M_cell, elev)

            # Surface moisture
            evap = compute_ocean_evaporation_rate(T_cell, elev)

            elev_meters = max(0.0, elev) * ELEVATION_SCALE
            T_at_altitude = T_cell - LAPSE_RATE * elev_meters

            drying_factor = compute_descent_saturation_multiplier(descent)
            M_sat_effective = compute_saturation_moisture_at_temperature(T_at_altitude) * drying_factor

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

    # Use Tsit5 (explicit) with increased maxiters for long simulations
    # Implicit solvers are too slow for this problem size (Jacobian is 3200x3200)
    # The numerical fixes (positivity constraints, tuned parameters) help stability
    if callback === nothing
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1800.0, maxiters=10_000_000)
    else
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1800.0, maxiters=10_000_000, callback=callback)
    end

    return sol
end

"""
    run_simulation_with_full_twolayer_atmosphere(moon::MoonBody2D, hours, T0, M0; U0=nothing, M_up0=nothing, T_up0=nothing, callback=nothing)

Run a full two-layer atmosphere simulation with upper temperature dynamics.
"""
function run_simulation_with_full_twolayer_atmosphere(moon::MoonBody2D, hours::Real,
                                                       T0::AbstractMatrix{<:Real},
                                                       M0::AbstractMatrix{<:Real};
                                                       U0::Union{AbstractMatrix{<:Real},Nothing}=nothing,
                                                       M_up0::Union{AbstractMatrix{<:Real},Nothing}=nothing,
                                                       T_up0::Union{AbstractMatrix{<:Real},Nothing}=nothing,
                                                       callback=nothing)
    n_lat, n_lon = moon.n_lat, moon.n_lon

    if U0 === nothing
        U0 = fill(U_INITIAL, n_lat, n_lon)
    end
    if M_up0 === nothing
        M_up0 = fill(M_UP_INITIAL, n_lat, n_lon)
    end
    if T_up0 === nothing
        T_up0 = fill(T_UP_INITIAL, n_lat, n_lon)
    end

    u0 = pack_full_twolayer_state(T0, M0, U0, M_up0, T_up0)
    tspan = (0.0, hours * 3600.0)

    prob = ODEProblem(derivatives_full_twolayer_atmosphere!, u0, tspan, moon)

    if callback === nothing
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1800.0, maxiters=10_000_000)
    else
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1800.0, maxiters=10_000_000, callback=callback)
    end

    return sol
end
