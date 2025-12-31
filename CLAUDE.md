# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Hot Moon Climate Simulator - a Julia-based scientific simulation for modeling 2D/1D climate dynamics on an exoplanet moon. Uses energy balance ODEs (not full fluid dynamics) to simulate temperature, moisture, and atmospheric circulation patterns.

## Build & Run Commands

```bash
# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run simulations
julia --project=. scripts/run_1d.jl              # 1D latitude-only
julia --project=. scripts/run_2d.jl              # 2D temperature
julia --project=. scripts/run_2d_moisture.jl     # 2D with moisture coupling
julia --project=. scripts/run_2d_twolayer.jl     # 2D with two-layer atmosphere

# With multithreading (recommended for 2D)
julia --project=. -t 4 scripts/run_2d.jl
```

## Architecture

```
HotMoon.jl (main module)
├── types.jl       → AbstractMoonBody, MoonBody1D, MoonBody2D
├── constants.jl   → Physical parameters (orbital, radiative, transport)
├── geometry.jl    → Solar geometry, zenith angles, eclipse detection
├── physics.jl     → Core physics: albedo, greenhouse, moisture, biomes
├── solver_1d.jl   → 1D ODE system (latitude bands)
├── solver_2d.jl   → 2D ODE system + moisture/two-layer variants
├── visualization.jl → Plotting, animation, analysis
└── performance.jl → Threading utilities
```

**Design Pattern:** Multiple dispatch on moon type enables unified `run_simulation(moon, hours, T0)` that auto-selects correct solver.

## Key Types

- `MoonBody1D(n_lat)` - Quick simulations with latitude bands only
- `MoonBody2D(n_lat, n_lon; seed, sea_level)` - Full 2D grid with terrain
- Unified constructor: `HotMoonBody(18)` → 1D, `HotMoonBody(90, 180)` → 2D

## Simulation Variants

1. **Temperature-only**: State = T[i,j]
2. **Temperature + Moisture**: State = [T; M] (evaporation-precipitation)
3. **Two-Layer Atmosphere**: State = [T; M; U; M_up] (vertical circulation)

## Important Patterns

**Smooth transitions (tanh)**: All feedbacks use tanh for numerical stability - ice-albedo, convection activation, biome blending.

**State vector flattening**: ODEs work with flat vectors. 2D indexing: `idx = (i-1) * n_lon + j`. Moisture/two-layer concatenate systems: `[T; M; U; M_up]`.

**Longitude wrapping**: Always use `mod1(j ± 1, n_lon)` for periodic boundary.

**Zenith returning nothing**: When sun below horizon, `get_zenith` returns `nothing` (not NaN). Check with `if zenith === nothing`.

**Pre-computed transport**: Terrain-based transport coefficients calculated at initialization and cached in moon structure.

**Time units**: Internal = seconds; user-facing = hours; save interval = 30 minutes.

## Physics Notes

- `ROTATION_PERIOD = 28 hrs`
- `SOLAR_CONSTANT = 2000 W/m²` (1.47× Earth)
- Heat transport: conduction (15 W/m²·K⁻¹) + convection (30 W/m²·K⁻¹ above 270K)
- Moisture barriers: Mountains block moisture (MOISTURE_BARRIER_STRENGTH = 8.0) more than heat (4.0)
- Elevation: Normalized -1 to 1 (ocean to mountains)
