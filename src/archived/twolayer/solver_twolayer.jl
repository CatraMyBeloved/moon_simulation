"""
Two-Layer Atmosphere Solvers (ARCHIVED)

This file contains the two-layer atmosphere ODE solvers that have been archived
due to parameter tuning challenges. The implementation includes:

1. Basic two-layer: T, M, U (upper mass), M_up (upper moisture)
2. Full two-layer: adds T_up (upper temperature dynamics)

To restore: include this file after loading the main HotMoon module,
and also include constants_twolayer.jl and atmosphere.jl from this folder.

Known issues:
- Upper layer precipitation caused runaway cooling
- Evaporation threshold (290K) created positive feedback loops
- Parameter sensitivity made equilibrium hard to achieve
"""

using DifferentialEquations

# =============================================================================
# Two-Layer Atmosphere Solver (4 fields: T, M, U, M_up)
# =============================================================================

"""
    derivatives_twolayer_atmosphere!(du, u, moon::MoonBody2D, t)

Two-layer atmosphere ODE system.
State vector layout: [T; M; U; M_up] where each is flattened n_lat × n_lon.

Physical processes:
- Surface: solar heating, IR cooling, evaporation, precipitation, transport
- Upper layer: ascent from surface, descent, poleward transport
- Moisture: lifted during ascent, precipitates from upper layer
"""
function derivatives_twolayer_atmosphere!(du, u, moon::MoonBody2D, t)
    n_lat, n_lon = moon.n_lat, moon.n_lon
    n_cells = n_lat * n_lon

    # Unpack state
    state = unpack_twolayer_state(u, moon)
    derivs = unpack_twolayer_derivatives(du, moon)

    T, M, U, M_up = state.T, state.M, state.U, state.M_up
    dT, dM, dU, dM_up = derivs.dT, derivs.dM, derivs.dU, derivs.dM_up

    # Get caches
    heat_transport = reshape(moon._transport_cache, n_lat, n_lon)
    moist_transport = reshape(moon._moisture_cache, n_lat, n_lon)
    upper_transport = reshape(moon._upper_transport_cache, n_lat, n_lon)
    upper_moist_transport = reshape(moon._upper_moist_cache, n_lat, n_lon)

    fill!(heat_transport, 0.0)
    fill!(moist_transport, 0.0)
    fill!(upper_transport, 0.0)
    fill!(upper_moist_transport, 0.0)

    # Surface transport
    foreach_cell(n_lat, n_lon) do i, j
        @inbounds compute_cell_heat_and_moisture_transport!(heat_transport, moist_transport, T, M, moon, i, j)
    end

    # Upper layer transport (mass and moisture)
    foreach_cell(n_lat, n_lon) do i, j
        @inbounds compute_upper_layer_transport!(upper_transport, upper_moist_transport, U, M_up, moon, i, j)
    end

    # Main physics loop
    foreach_cell(n_lat, n_lon) do i, j
        @inbounds begin
            T_cell = T[i, j]
            M_cell = max(0.0, M[i, j])
            U_cell = max(U_FLOOR, U[i, j])
            M_up_cell = max(0.0, M_up[i, j])
            lat = moon.latitudes[i]
            lon = moon.longitudes[j]
            elev = moon.elevation[i, j]

            # Solar input
            Q_in = get_solar_2d(t, lat, lon, T_cell)

            # Greenhouse effect
            τ_IR = compute_total_ir_optical_depth(M_cell)
            τ_effective = compute_effective_optical_depth_at_elevation(τ_IR, elev)
            transmissivity = compute_ir_transmissivity_from_optical_depth(τ_effective)
            Q_out = compute_outgoing_longwave_radiation(T_cell, transmissivity)

            # Heat capacity
            heat_cap = compute_biome_blended_heat_capacity(T_cell, M_cell, elev)

            # Surface moisture processes
            evap = compute_ocean_evaporation_rate(T_cell, elev)
            precip = compute_orographic_precipitation_rate(M_cell, T_cell, elev)

            # Vertical circulation
            ascent = compute_ascent_rate(T_cell, M_cell, lat)
            descent = compute_descent_rate(U_cell, lat)

            # Moisture exchange
            moisture_lift = compute_moisture_lift(M_cell, T_cell, ascent)
            moisture_fall = compute_moisture_fall(M_up_cell, descent)

            # Upper layer precipitation
            T_up_approx = T_cell - UPPER_LAYER_TEMP_DROP
            precip_upper = compute_upper_precipitation(M_up_cell, T_up_approx, U_cell)

            # Latent heat
            Q_latent = compute_latent_heat_flux(precip + precip_upper, evap)

            # Floor restoration for U
            u_restoration = U_cell < U_FLOOR ? U_FLOOR_RESTORATION_RATE * (U_FLOOR - U_cell) : 0.0

            # Derivatives
            dT[i, j] = (Q_in - Q_out + heat_transport[i, j] + Q_latent) / heat_cap
            dM[i, j] = evap - precip + moisture_fall + precip_upper - moisture_lift + moist_transport[i, j]
            dU[i, j] = ascent - descent + upper_transport[i, j] + u_restoration
            dM_up[i, j] = moisture_lift - moisture_fall - precip_upper + upper_moist_transport[i, j]
        end
    end
end

# =============================================================================
# Full Two-Layer Atmosphere Solver (5 fields: T, M, U, M_up, T_up)
# =============================================================================

"""
    derivatives_full_twolayer_atmosphere!(du, u, moon::MoonBody2D, t)

Full two-layer atmosphere with upper temperature dynamics.
State vector layout: [T; M; U; M_up; T_up]
"""
function derivatives_full_twolayer_atmosphere!(du, u, moon::MoonBody2D, t)
    n_lat, n_lon = moon.n_lat, moon.n_lon
    n_cells = n_lat * n_lon

    # Unpack state
    state = unpack_full_twolayer_state(u, moon)
    derivs = unpack_full_twolayer_derivatives(du, moon)

    T, M, U, M_up, T_up = state.T, state.M, state.U, state.M_up, state.T_up
    dT, dM, dU, dM_up, dT_up = derivs.dT, derivs.dM, derivs.dU, derivs.dM_up, derivs.dT_up

    # Get caches
    heat_transport = reshape(moon._transport_cache, n_lat, n_lon)
    moist_transport = reshape(moon._moisture_cache, n_lat, n_lon)
    upper_transport = reshape(moon._upper_transport_cache, n_lat, n_lon)
    upper_moist_transport = reshape(moon._upper_moist_cache, n_lat, n_lon)
    upper_heat_transport = reshape(moon._upper_heat_cache, n_lat, n_lon)

    fill!(heat_transport, 0.0)
    fill!(moist_transport, 0.0)
    fill!(upper_transport, 0.0)
    fill!(upper_moist_transport, 0.0)
    fill!(upper_heat_transport, 0.0)

    # Surface transport
    foreach_cell(n_lat, n_lon) do i, j
        @inbounds compute_cell_heat_and_moisture_transport!(heat_transport, moist_transport, T, M, moon, i, j)
    end

    # Upper layer transport (mass, moisture, and heat)
    foreach_cell(n_lat, n_lon) do i, j
        @inbounds compute_upper_layer_full_transport!(upper_transport, upper_moist_transport, upper_heat_transport,
                                                       U, M_up, T_up, moon, i, j)
    end

    # Main physics loop
    foreach_cell(n_lat, n_lon) do i, j
        @inbounds begin
            T_cell = T[i, j]
            M_cell = max(0.0, M[i, j])
            U_cell = max(U_FLOOR, U[i, j])
            M_up_cell = max(0.0, M_up[i, j])
            T_up_cell = clamp(T_up[i, j], T_UP_FLOOR, T_UP_CEILING)
            lat = moon.latitudes[i]
            lon = moon.longitudes[j]
            elev = moon.elevation[i, j]

            # Solar input
            Q_in = get_solar_2d(t, lat, lon, T_cell)

            # Greenhouse effect
            τ_IR = compute_total_ir_optical_depth(M_cell)
            τ_effective = compute_effective_optical_depth_at_elevation(τ_IR, elev)
            transmissivity = compute_ir_transmissivity_from_optical_depth(τ_effective)
            Q_out = compute_outgoing_longwave_radiation(T_cell, transmissivity)

            # Heat capacity
            heat_cap = compute_biome_blended_heat_capacity(T_cell, M_cell, elev)

            # Surface moisture processes
            evap = compute_ocean_evaporation_rate(T_cell, elev)
            precip = compute_orographic_precipitation_rate(M_cell, T_cell, elev)

            # Vertical circulation with temperature feedback
            ascent = compute_ascent_rate_with_temperature(T_cell, M_cell, T_up_cell, lat)
            descent = compute_descent_rate_with_temperature(U_cell, T_up_cell, lat)

            # Moisture exchange
            moisture_lift = compute_moisture_lift(M_cell, T_cell, ascent)
            moisture_fall = compute_moisture_fall(M_up_cell, descent)

            # Upper layer precipitation
            precip_upper = compute_upper_precipitation(M_up_cell, T_up_cell, U_cell)

            # Latent heat
            Q_latent = compute_latent_heat_flux(precip + precip_upper, evap)

            # Floor restoration for U
            u_restoration = U_cell < U_FLOOR ? U_FLOOR_RESTORATION_RATE * (U_FLOOR - U_cell) : 0.0

            # Upper layer temperature dynamics
            heat_lift = compute_heat_lift(T_cell, T_up_cell, ascent)
            heat_fall = compute_heat_fall(T_up_cell, descent)
            radiative_cooling = compute_upper_radiative_cooling(T_up_cell)

            # Derivatives
            dT[i, j] = (Q_in - Q_out + heat_transport[i, j] + Q_latent) / heat_cap
            dM[i, j] = evap - precip + moisture_fall + precip_upper - moisture_lift + moist_transport[i, j]
            dU[i, j] = ascent - descent + upper_transport[i, j] + u_restoration
            dM_up[i, j] = moisture_lift - moisture_fall - precip_upper + upper_moist_transport[i, j]
            dT_up[i, j] = (heat_lift - heat_fall + upper_heat_transport[i, j] + radiative_cooling) / UPPER_LAYER_HEAT_CAPACITY
        end
    end
end

# =============================================================================
# Upper Layer Transport Helpers
# =============================================================================

"""
    compute_upper_layer_transport!(upper_transport, upper_moist_transport, U, M_up, moon, i, j)

Compute upper layer mass and moisture transport for basic two-layer model.
"""
@inline function compute_upper_layer_transport!(upper_transport, upper_moist_transport, U, M_up, moon, i, j)
    n_lat, n_lon = moon.n_lat, moon.n_lon
    lat = moon.latitudes[i]

    U_self = max(U_FLOOR, U[i, j])
    M_up_self = max(0.0, M_up[i, j])

    # Meridional transport (poleward)
    if i > 1
        U_nb = max(U_FLOOR, U[i-1, j])
        M_up_nb = max(0.0, M_up[i-1, j])
        flow = UPPER_MERIDIONAL_COEFF * (U_nb - U_self)
        upper_transport[i, j] += flow
        upper_moist_transport[i, j] += flow > 0 ? flow * M_up_nb / U_nb : flow * M_up_self / U_self
    end
    if i < n_lat
        U_nb = max(U_FLOOR, U[i+1, j])
        M_up_nb = max(0.0, M_up[i+1, j])
        flow = UPPER_MERIDIONAL_COEFF * (U_nb - U_self)
        upper_transport[i, j] += flow
        upper_moist_transport[i, j] += flow > 0 ? flow * M_up_nb / U_nb : flow * M_up_self / U_self
    end

    # Zonal transport
    j_east = mod1(j + 1, n_lon)
    j_west = mod1(j - 1, n_lon)

    U_east = max(U_FLOOR, U[i, j_east])
    U_west = max(U_FLOOR, U[i, j_west])
    M_up_east = max(0.0, M_up[i, j_east])
    M_up_west = max(0.0, M_up[i, j_west])

    # Pressure-driven zonal flow
    flow_east = UPPER_ZONAL_COEFF * (U_east - U_self)
    flow_west = UPPER_ZONAL_COEFF * (U_west - U_self)

    # Westerly bias (Coriolis-like)
    westerly = WESTERLY_BIAS_STRENGTH * cos(deg2rad(lat))
    flow_east += westerly
    flow_west -= westerly

    upper_transport[i, j] += flow_east + flow_west
    upper_moist_transport[i, j] += (flow_east > 0 ? flow_east * M_up_east / U_east : flow_east * M_up_self / U_self)
    upper_moist_transport[i, j] += (flow_west > 0 ? flow_west * M_up_west / U_west : flow_west * M_up_self / U_self)
end

"""
    compute_upper_layer_full_transport!(upper_transport, upper_moist_transport, upper_heat_transport, U, M_up, T_up, moon, i, j)

Compute upper layer mass, moisture, and heat transport for full two-layer model.
"""
@inline function compute_upper_layer_full_transport!(upper_transport, upper_moist_transport, upper_heat_transport,
                                                      U, M_up, T_up, moon, i, j)
    n_lat, n_lon = moon.n_lat, moon.n_lon
    lat = moon.latitudes[i]

    U_self = max(U_FLOOR, U[i, j])
    M_up_self = max(0.0, M_up[i, j])
    T_up_self = T_up[i, j]

    # Meridional transport
    if i > 1
        U_nb = max(U_FLOOR, U[i-1, j])
        M_up_nb = max(0.0, M_up[i-1, j])
        T_up_nb = T_up[i-1, j]
        flow = UPPER_MERIDIONAL_COEFF * (U_nb - U_self)
        upper_transport[i, j] += flow
        upper_moist_transport[i, j] += flow > 0 ? flow * M_up_nb / U_nb : flow * M_up_self / U_self
        # Temperature is intensive - don't divide by U
        upper_heat_transport[i, j] += flow > 0 ? flow * T_up_nb : flow * T_up_self
    end
    if i < n_lat
        U_nb = max(U_FLOOR, U[i+1, j])
        M_up_nb = max(0.0, M_up[i+1, j])
        T_up_nb = T_up[i+1, j]
        flow = UPPER_MERIDIONAL_COEFF * (U_nb - U_self)
        upper_transport[i, j] += flow
        upper_moist_transport[i, j] += flow > 0 ? flow * M_up_nb / U_nb : flow * M_up_self / U_self
        upper_heat_transport[i, j] += flow > 0 ? flow * T_up_nb : flow * T_up_self
    end

    # Zonal transport
    j_east = mod1(j + 1, n_lon)
    j_west = mod1(j - 1, n_lon)

    U_east = max(U_FLOOR, U[i, j_east])
    U_west = max(U_FLOOR, U[i, j_west])
    M_up_east = max(0.0, M_up[i, j_east])
    M_up_west = max(0.0, M_up[i, j_west])
    T_up_east = T_up[i, j_east]
    T_up_west = T_up[i, j_west]

    flow_east = UPPER_ZONAL_COEFF * (U_east - U_self)
    flow_west = UPPER_ZONAL_COEFF * (U_west - U_self)

    westerly = WESTERLY_BIAS_STRENGTH * cos(deg2rad(lat))
    flow_east += westerly
    flow_west -= westerly

    upper_transport[i, j] += flow_east + flow_west
    upper_moist_transport[i, j] += (flow_east > 0 ? flow_east * M_up_east / U_east : flow_east * M_up_self / U_self)
    upper_moist_transport[i, j] += (flow_west > 0 ? flow_west * M_up_west / U_west : flow_west * M_up_self / U_self)
    upper_heat_transport[i, j] += (flow_east > 0 ? flow_east * T_up_east : flow_east * T_up_self)
    upper_heat_transport[i, j] += (flow_west > 0 ? flow_west * T_up_west : flow_west * T_up_self)
end

# =============================================================================
# Simulation Runners
# =============================================================================

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
