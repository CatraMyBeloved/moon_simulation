"""
1D ODE solver for latitude-band climate model
"""

"""
    derivatives_1d!(dT, temps, moon::MoonBody1D, t)

Calculate temperature derivatives for each latitude band.

This is the core ODE system: dT/dt = (heating - cooling + transport) / C
"""
function derivatives_1d!(dT, temps, moon::MoonBody1D, t)
    n = moon.n_lat
    transport = moon._transport_cache

    fill!(transport, 0.0)
    for i in 2:n
        # Non-linear transport: use average temperature to determine coefficient
        T_avg = 0.5 * (temps[i-1] + temps[i])
        k = get_transport_coefficient(T_avg)
        flow = k * (temps[i-1] - temps[i])
        transport[i] += flow
        transport[i-1] -= flow
    end

    for i in 1:n
        Q_in = get_solar(t, moon.latitudes[i], temps[i])

        greenhouse = get_greenhouse(temps[i])

        Q_out = EMISSIVITY * STEFAN_BOLTZMANN * temps[i]^4 * (1 - greenhouse)

        # Use latitude-dependent heat capacity
        heat_cap = get_heat_capacity(moon.latitudes[i])
        dT[i] = (Q_in - Q_out + transport[i]) / heat_cap
    end
end

"""
    run_simulation(moon::MoonBody1D, hours, T0)

Run a 1D climate simulation.

# Arguments
- `moon`: MoonBody1D structure
- `hours`: Simulation duration in hours
- `T0`: Initial temperatures (Kelvin) for each band

# Returns
- Solution object from DifferentialEquations.jl
"""
function run_simulation(moon::MoonBody1D, hours::Real, T0::AbstractVector{<:Real})
    tspan = (0.0, hours * 3600)

    prob = ODEProblem(derivatives_1d!, T0, tspan, moon)

    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8, saveat=1800.0)

    return sol
end
