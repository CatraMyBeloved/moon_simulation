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

        # IR optical depth and transmissivity (dry atmosphere for 1D)
        τ_IR = get_ir_optical_depth()
        transmissivity = get_ir_transmissivity(τ_IR)

        Q_out = EMISSIVITY * STEFAN_BOLTZMANN * temps[i]^4 * transmissivity

        # Use latitude-dependent heat capacity
        heat_cap = get_heat_capacity(moon.latitudes[i])
        dT[i] = (Q_in - Q_out + transport[i]) / heat_cap
    end
end

"""
    run_simulation(moon::MoonBody1D, hours, T0; callback=nothing)

Run a 1D climate simulation.

# Arguments
- `moon`: MoonBody1D structure
- `hours`: Simulation duration in hours
- `T0`: Initial temperatures (Kelvin) for each band
- `callback`: Optional callback for progress reporting (default: nothing)

# Returns
- Solution object from DifferentialEquations.jl
"""
function run_simulation(moon::MoonBody1D, hours::Real, T0::AbstractVector{<:Real}; callback=nothing)
    tspan = (0.0, hours * 3600)

    prob = ODEProblem(derivatives_1d!, T0, tspan, moon)

    if callback === nothing
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8, saveat=1800.0)
    else
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8, saveat=1800.0, callback=callback)
    end

    return sol
end

"""
    make_progress_callback(total_hours; update_interval_hours=500)

Create a callback that prints simulation progress.

# Arguments
- `total_hours`: Total simulation duration in hours
- `update_interval_hours`: How often to print progress (default: 500 hours)

# Returns
- PeriodicCallback suitable for passing to run_simulation
"""
function make_progress_callback(total_hours::Real; update_interval_hours::Real=100)
    total_sec = total_hours * 3600
    interval_sec = update_interval_hours * 3600

    # Print initial message
    println("  [Progress updates every $(round(Int, update_interval_hours)) hours]")

    callback = PeriodicCallback(interval_sec) do integrator
        pct = 100 * integrator.t / total_sec
        hr = round(Int, integrator.t / 3600)
        println("  Progress: $(lpad(round(Int, pct), 3))% ($(hr) / $(round(Int, total_hours)) hours)")
    end

    return callback
end
