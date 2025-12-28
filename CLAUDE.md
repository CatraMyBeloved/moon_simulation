# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Hot Moon Climate Simulator - a 1D/2D climate simulation for a moon orbiting a gas giant in a red giant star system. Uses energy balance ODEs to model temperature patterns, day/night cycles, eclipses, and topographic effects without fluid dynamics.

**Setting**: The moon is *not* tidally locked - it rotates independently (28-hour day) while orbiting its gas giant (64-hour orbit). This creates eclipse events that drift through local time, hitting different hours each orbit. The incommensurable periods (GCD=4, LCM=448) produce a 16-day/7-orbit cycle before patterns repeat.

**Purpose**: Foundation for worldbuilding - the simulation informs climate zones, habitable regions, and natural calendar systems for a fantasy world with its own cultures, timekeeping, and inhabitants.

## Commands

```bash
# Run 1D simulation (latitude bands only)
julia --project=. scripts/run_1d.jl

# Run 2D simulation (latitude × longitude grid)
julia --project=. scripts/run_2d.jl

# Run tests
julia --project=. tests/test_physics.jl

# Install dependencies (uses Project.toml)
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Architecture

The codebase uses a Julia module (`HotMoon`) defined in `src/HotMoon.jl`. Scripts load it via:
```julia
include(joinpath(@__DIR__, "../src/HotMoon.jl"))
using .HotMoon
```

### Multiple Dispatch Pattern

The codebase uses Julia's multiple dispatch for 1D/2D support:

```julia
# Type hierarchy
abstract type AbstractMoonBody end
struct MoonBody1D <: AbstractMoonBody  # latitude bands only
struct MoonBody2D <: AbstractMoonBody  # latitude × longitude grid

# Unified constructor
moon1d = HotMoonBody(18)        # Returns MoonBody1D with 18 latitude bands
moon2d = HotMoonBody(18, 36)    # Returns MoonBody2D with 18×36 grid

# Unified simulation function (dispatches on type)
sol = run_simulation(moon, hours, T0)  # Works for both 1D and 2D
```

### Core Physics Pipeline

1. **types.jl** - `AbstractMoonBody` with `MoonBody1D` and `MoonBody2D` implementations
2. **constants.jl** - Physical parameters including zone-based heat capacities and non-linear transport
3. **physics.jl** - Temperature-dependent albedo, greenhouse effect, transport coefficients, and terrain generation
4. **geometry.jl** - Solar zenith angles (1D and 2D), eclipse detection
5. **solver_1d.jl** - 1D ODE system: `dT/dt = (Q_in - Q_out + transport) / C`
6. **solver_2d.jl** - 2D ODE system with 4-direction transport and longitude wrapping
7. **visualization.jl** - Plots for both 1D and 2D simulations (global mean, profiles, heatmaps, Hovmoeller diagrams)

### Key Physical Model Details

- **Energy balance**: Solar heating vs Stefan-Boltzmann emission modulated by greenhouse factor
- **Heat transport**: Non-linear diffusion with convection activation (temperature-dependent)
- **Ice-albedo feedback**: Higher albedo below freezing (tanh transition around 273K)
- **Water vapor feedback**: Stronger greenhouse effect at higher temperatures
- **Zone-based thermal inertia**: Equatorial zones have 10x higher heat capacity than polar
- **Terrain effects**: Fractal noise elevation affects heat capacity and transport (mountains block, oceans enhance)

### Timing Parameters

- Rotation period: 28 hours
- Orbital period: 64 hours
- Eclipse duration: 4 hours (occurs at orbital_phase = 0.5)

## Planned but Not Implemented

- Longitude-dependent heat capacity for oceans/continents
- Pluto notebooks for interactive exploration
