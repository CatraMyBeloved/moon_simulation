"""
1D ODE solver for latitude-band climate model.
"""

using DifferentialEquations

"""
    derivatives_1d!(dT, temps, moon::MoonBody1D, t)

Calculate temperature derivatives for each latitude band.
This is the core ODE system: dT/dt = (heating - cooling + transport) / C
"""
function derivatives_1d!(dT, temps, moon::MoonBody1D, t)
    n = moon.n_lat
    transport = moon._transport_cache

    # Transport loop - has sequential dependencies, cannot parallelize
    fill!(transport, 0.0)
    @inbounds for i in 2:n
        T_avg = 0.5 * (temps[i-1] + temps[i])
        k = compute_heat_transport_coefficient(T_avg)
        flow = k * (temps[i-1] - temps[i])
        transport[i] += flow
        transport[i-1] -= flow
    end

    # Physics loop - embarrassingly parallel
    if should_thread(n)
        @inbounds Threads.@threads for i in 1:n
            Q_in = get_solar(t, moon.latitudes[i], temps[i])

            τ_IR = compute_total_ir_optical_depth()
            transmissivity = compute_ir_transmissivity_from_optical_depth(τ_IR)

            Q_out = compute_outgoing_longwave_radiation(temps[i], transmissivity)

            heat_cap = compute_latitude_based_heat_capacity(moon.latitudes[i])
            dT[i] = (Q_in - Q_out + transport[i]) / heat_cap
        end
    else
        @inbounds @simd for i in 1:n
            Q_in = get_solar(t, moon.latitudes[i], temps[i])

            τ_IR = compute_total_ir_optical_depth()
            transmissivity = compute_ir_transmissivity_from_optical_depth(τ_IR)

            Q_out = compute_outgoing_longwave_radiation(temps[i], transmissivity)

            heat_cap = compute_latitude_based_heat_capacity(moon.latitudes[i])
            dT[i] = (Q_in - Q_out + transport[i]) / heat_cap
        end
    end
end

"""
    run_simulation(moon::MoonBody1D, hours, T0; callback=nothing)

Run a 1D climate simulation.
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
    make_progress_callback(total_hours; update_interval_hours=100)

Create a callback that prints simulation progress.
"""
function make_progress_callback(total_hours::Real; update_interval_hours::Real=100)
    total_sec = total_hours * 3600
    interval_sec = update_interval_hours * 3600

    println("  [Progress updates every $(round(Int, update_interval_hours)) hours]")

    callback = PeriodicCallback(interval_sec) do integrator
        pct = 100 * integrator.t / total_sec
        hr = round(Int, integrator.t / 3600)
        println("  Progress: $(lpad(round(Int, pct), 3))% ($(hr) / $(round(Int, total_hours)) hours)")
    end

    return callback
end
