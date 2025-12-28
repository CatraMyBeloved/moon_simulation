"""
2D ODE solver for latitude-longitude climate model
"""

"""
    derivatives_2d!(dT, temps, moon::MoonBody2D, t)

Calculate temperature derivatives for each grid cell.

The ODE solver passes flat vectors, so we reshape to 2D for the calculation.
For each cell: dT/dt = (Q_solar - Q_radiated + Q_transport) / heat_capacity
"""
function derivatives_2d!(dT, temps, moon::MoonBody2D, t)
    n_lat = moon.n_lat
    n_lon = moon.n_lon

    T = reshape(temps, n_lat, n_lon)
    dT_2d = reshape(dT, n_lat, n_lon)
    transport = reshape(moon._transport_cache, n_lat, n_lon)


    fill!(transport, 0.0)

    for i in 1:n_lat
        for j in 1:n_lon
            T_self = T[i, j]

            # North neighbor
            if i > 1
                T_neighbor = T[i-1, j]
                k = get_transport_coefficient(0.5 * (T_self + T_neighbor))
                k_dir = k * moon.transport_coeffs[i, j, 1]
                transport[i, j] += k_dir * (T_neighbor - T_self)
            end

            # South neighbor
            if i < n_lat
                T_neighbor = T[i+1, j]
                k = get_transport_coefficient(0.5 * (T_self + T_neighbor))
                k_dir = k * moon.transport_coeffs[i, j, 2]
                transport[i, j] += k_dir * (T_neighbor - T_self)
            end

            # East neighbor (wraps)
            j_east = mod1(j + 1, n_lon)
            T_neighbor = T[i, j_east]
            k = get_transport_coefficient(0.5 * (T_self + T_neighbor))
            k_dir = k * moon.transport_coeffs[i, j, 3]
            transport[i, j] += k_dir * (T_neighbor - T_self)

            # West neighbor (wraps)
            j_west = mod1(j - 1, n_lon)
            T_neighbor = T[i, j_west]
            k = get_transport_coefficient(0.5 * (T_self + T_neighbor))
            k_dir = k * moon.transport_coeffs[i, j, 4]
            transport[i, j] += k_dir * (T_neighbor - T_self)
        end
    end

    for i in 1:n_lat
        for j in 1:n_lon
            T_cell = T[i, j]
            lat = moon.latitudes[i]
            lon = moon.longitudes[j]

            # Solar heating: depends on time, position, and albedo (temperature-dependent)
            Q_in = get_solar_2d(t, lat, lon, T_cell)

            greenhouse = get_greenhouse(T_cell)
            Q_out = EMISSIVITY * STEFAN_BOLTZMANN * T_cell^4 * (1 - greenhouse)

            # Use elevation to determine land vs ocean heat capacity
            elev = moon.elevation[i, j]
            heat_cap = get_heat_capacity(lat, lon, elev)

            dT_2d[i, j] = (Q_in - Q_out + transport[i, j]) / heat_cap
        end
    end
end

"""
    run_simulation(moon::MoonBody2D, hours, T0)

Run a 2D climate simulation.

# Arguments
- `moon`: MoonBody2D structure
- `hours`: Simulation duration in hours
- `T0`: Initial temperatures (Kelvin) as n_lat × n_lon matrix

# Returns
- Solution object from DifferentialEquations.jl (states are flattened)
"""
function run_simulation(moon::MoonBody2D, hours::Real, T0::AbstractMatrix{<:Real})
    tspan = (0.0, hours * 3600)

    T0_flat = vec(T0)

    prob = ODEProblem(derivatives_2d!, T0_flat, tspan, moon)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8, saveat=1800.0)

    return sol
end

# =============================================================================
# Coupled Temperature-Moisture Solver
# =============================================================================

"""
    derivatives_2d_moisture!(du, u, moon::MoonBody2D, t)

Coupled temperature-moisture ODE system.
State vector layout: [T₁₁, ..., Tₙₘ, M₁₁, ..., Mₙₘ]

Temperature equation: dT/dt = (Q_solar - Q_radiated + Q_transport) / heat_capacity
Moisture equation: dM/dt = evaporation - precipitation + moisture_transport
"""
function derivatives_2d_moisture!(du, u, moon::MoonBody2D, t)
    n_lat = moon.n_lat
    n_lon = moon.n_lon
    n_cells = n_lat * n_lon

    # Split state vector into temperature and moisture
    T = reshape(@view(u[1:n_cells]), n_lat, n_lon)
    M = reshape(@view(u[n_cells+1:end]), n_lat, n_lon)

    dT = reshape(@view(du[1:n_cells]), n_lat, n_lon)
    dM = reshape(@view(du[n_cells+1:end]), n_lat, n_lon)

    # Transport caches
    heat_transport = reshape(moon._transport_cache, n_lat, n_lon)
    moisture_transport = reshape(moon._moisture_cache, n_lat, n_lon)

    fill!(heat_transport, 0.0)
    fill!(moisture_transport, 0.0)

    # ========== TRANSPORT LOOPS ==========

    for i in 1:n_lat
        for j in 1:n_lon
            T_self = T[i, j]
            M_self = M[i, j]

            # North neighbor
            if i > 1
                T_nb, M_nb = T[i-1, j], M[i-1, j]

                # Heat transport
                k_heat = get_transport_coefficient(0.5 * (T_self + T_nb))
                k_heat_dir = k_heat * moon.transport_coeffs[i, j, 1]
                heat_transport[i, j] += k_heat_dir * (T_nb - T_self)

                # Moisture transport
                k_moist_dir = MOISTURE_DIFFUSION * moon.moisture_transport_coeffs[i, j, 1]
                moisture_transport[i, j] += k_moist_dir * (M_nb - M_self)
            end

            # South neighbor
            if i < n_lat
                T_nb, M_nb = T[i+1, j], M[i+1, j]

                k_heat = get_transport_coefficient(0.5 * (T_self + T_nb))
                k_heat_dir = k_heat * moon.transport_coeffs[i, j, 2]
                heat_transport[i, j] += k_heat_dir * (T_nb - T_self)

                k_moist_dir = MOISTURE_DIFFUSION * moon.moisture_transport_coeffs[i, j, 2]
                moisture_transport[i, j] += k_moist_dir * (M_nb - M_self)
            end

            # East neighbor (wraps)
            j_east = mod1(j + 1, n_lon)
            T_nb, M_nb = T[i, j_east], M[i, j_east]

            k_heat = get_transport_coefficient(0.5 * (T_self + T_nb))
            k_heat_dir = k_heat * moon.transport_coeffs[i, j, 3]
            heat_transport[i, j] += k_heat_dir * (T_nb - T_self)

            k_moist_dir = MOISTURE_DIFFUSION * moon.moisture_transport_coeffs[i, j, 3]
            moisture_transport[i, j] += k_moist_dir * (M_nb - M_self)

            # West neighbor (wraps)
            j_west = mod1(j - 1, n_lon)
            T_nb, M_nb = T[i, j_west], M[i, j_west]

            k_heat = get_transport_coefficient(0.5 * (T_self + T_nb))
            k_heat_dir = k_heat * moon.transport_coeffs[i, j, 4]
            heat_transport[i, j] += k_heat_dir * (T_nb - T_self)

            k_moist_dir = MOISTURE_DIFFUSION * moon.moisture_transport_coeffs[i, j, 4]
            moisture_transport[i, j] += k_moist_dir * (M_nb - M_self)
        end
    end

    # ========== BALANCE EQUATIONS ==========

    for i in 1:n_lat
        for j in 1:n_lon
            T_cell = T[i, j]
            M_cell = max(0.0, M[i, j])  # Ensure non-negative moisture
            lat = moon.latitudes[i]
            lon = moon.longitudes[j]
            elev = moon.elevation[i, j]

            # --- Temperature equation ---
            Q_in = get_solar_2d(t, lat, lon, T_cell)
            greenhouse = get_greenhouse(T_cell)
            Q_out = EMISSIVITY * STEFAN_BOLTZMANN * T_cell^4 * (1.0 - greenhouse)
            heat_cap = get_heat_capacity(lat, lon, elev)

            dT[i, j] = (Q_in - Q_out + heat_transport[i, j]) / heat_cap

            # --- Moisture equation ---
            evap = get_evaporation(T_cell, elev)
            precip = get_precipitation(M_cell, T_cell, elev)

            dM[i, j] = evap - precip + moisture_transport[i, j]
        end
    end
end

"""
    run_simulation_moisture(moon::MoonBody2D, hours, T0, M0)

Run a coupled temperature-moisture simulation.

# Arguments
- `moon`: MoonBody2D structure
- `hours`: Simulation duration in hours
- `T0`: Initial temperatures (Kelvin) as n_lat × n_lon matrix
- `M0`: Initial moisture (kg/m²) as n_lat × n_lon matrix

# Returns
- Solution object from DifferentialEquations.jl
  State layout: first n_lat*n_lon values are temperature, rest are moisture
"""
function run_simulation_moisture(moon::MoonBody2D, hours::Real,
                                  T0::AbstractMatrix{<:Real},
                                  M0::AbstractMatrix{<:Real})
    # Combine initial conditions: [T..., M...]
    u0 = vcat(vec(T0), vec(M0))

    tspan = (0.0, hours * 3600.0)

    prob = ODEProblem(derivatives_2d_moisture!, u0, tspan, moon)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8, saveat=1800.0)

    return sol
end
