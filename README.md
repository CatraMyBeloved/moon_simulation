# Hot Moon Climate Simulator

A scientific Julia framework for modeling 2D climate dynamics on an exoplanet moon. Simulates temperature, moisture, atmospheric circulation, and biome evolution using energy balance equations.

## Overview

Hot Moon simulates the climate physics of a moon orbiting a gas giant exoplanet:

- **28-hour rotation period** - slower than Earth, creating extended day/night cycles
- **64-hour orbital period** - moon completes an orbit around its planet in ~2.7 days
- **Periodic eclipses** - 4-hour eclipses when the moon passes behind its planet
- **Solar constant of 2000 W/m²** - 1.47× Earth's insolation
- **Procedural terrain** - fractal-generated mountains, oceans, and elevation-dependent effects

The simulator uses ordinary differential equations rather than full fluid dynamics, enabling efficient long-duration simulations (100,000+ hours) while capturing key climate feedbacks.

## Features

### Physics Engine
- **Radiative transfer** with ice-albedo feedback and water vapor greenhouse effect
- **Heat transport** via conduction and temperature-dependent convection
- **Moisture dynamics** with Clausius-Clapeyron saturation, evaporation, and orographic precipitation
- **Terrain effects** including rain shadows, lapse rate cooling, and directional transport barriers
- **Two-layer atmosphere** model simulating Hadley-like cells and subtropical dry zones

### Simulation Modes
| Mode | State Variables | Use Case |
|------|-----------------|----------|
| 1D Latitude | Temperature | Quick tests, educational |
| 2D Temperature | T[lat,lon] | Basic climate patterns |
| 2D + Moisture | T, M | Full water cycle |
| 2D Two-Layer | T, M, U, M_up | Atmospheric circulation |

### Biome System
12 distinct biome classes with unique thermal properties:

| Biome | Heat Capacity | Conditions |
|-------|---------------|------------|
| Ocean | 2.0×10⁷ J/m²/K | Elevation < 0 |
| Ice Sheet | 1.5×10⁷ | T < -10°C |
| Tropical Forest | 8.0×10⁶ | T > 20°C, wet |
| Hot Desert | 5.0×10⁵ | T > 20°C, dry |
| Mountain | 3.0×10⁵ | Elevation > 0.5 |

### Visualization
- Global temperature/moisture heatmaps
- Latitude-time Hovmöller diagrams
- Animated GIF evolution
- Biome distribution maps
- Precipitation patterns

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/moon_simulation.git
cd moon_simulation

# Install Julia dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

**Requirements:** Julia 1.9+ with packages: DifferentialEquations, Plots, JLD2

## Quick Start

### 1D Simulation (Latitude Bands)

```julia
using HotMoon

# Create 18-band latitude model
moon = HotMoonBody(18)
T0 = fill(285.0, 18)  # 12°C uniform start

# Run for 500 hours
sol = run_simulation(moon, 500, T0)

# Plot results
plot_1d_results(moon, sol)
```

### 2D Temperature Simulation

```julia
# Create 90×180 grid with terrain
moon = HotMoonBody(90, 180, seed=42, sea_level=-0.2)
T0 = fill(285.0, 90, 180)

# Run with 4 threads (recommended)
sol = run_simulation(moon, 5000, T0)

# Generate snapshot at final time
plot_2d_snapshot(moon, sol, sol.t[end])
```

### 2D with Moisture Coupling

```julia
moon = HotMoonBody(90, 180, seed=42)
T0, M0 = estimate_initial_temperature_and_moisture(moon)

sol = run_simulation_with_moisture(moon, 10000, T0, M0)

# Plot temperature, moisture, and precipitation
plot_climate_dashboard(moon, sol, sol.t[end])
```

### Two-Layer Atmosphere

```julia
moon = HotMoonBody(45, 90)
T0, M0 = estimate_initial_temperature_and_moisture(moon)
U0 = fill(1.0, 45, 90)      # Upper atmosphere mass
M_up0 = fill(0.05, 45, 90)  # Upper moisture

sol = run_simulation_with_twolayer_atmosphere(
    moon, 15000, T0, M0, U0=U0, M_up0=M_up0
)
```

## Command Line Usage

```bash
# Run simulations directly
julia --project=. scripts/run_1d.jl
julia --project=. scripts/run_2d.jl
julia --project=. scripts/run_2d_moisture.jl
julia --project=. scripts/run_2d_twolayer.jl

# Enable multithreading (recommended for 2D)
julia --project=. -t 4 scripts/run_2d_moisture.jl
```

## Project Structure

```
HotMoon/
├── src/
│   ├── HotMoon.jl          # Main module
│   ├── types.jl            # MoonBody1D, MoonBody2D
│   ├── constants.jl        # Physical parameters
│   ├── geometry.jl         # Solar angles, eclipses
│   ├── state.jl            # State vector utilities
│   │
│   ├── physics/
│   │   ├── radiative.jl    # Albedo, greenhouse, IR
│   │   ├── terrain.jl      # Fractal generation, transport
│   │   ├── moisture.jl     # Evaporation, precipitation
│   │   ├── biomes.jl       # Classification, heat capacity
│   │   └── atmosphere.jl   # Two-layer dynamics
│   │
│   ├── solvers/
│   │   ├── solver_1d.jl    # Latitude-band ODE
│   │   └── solver_2d.jl    # Full 2D variants
│   │
│   └── visualization/
│       ├── viz_1d.jl       # 1D plots
│       ├── viz_2d_*.jl     # 2D snapshots, timeseries
│       └── viz_animation.jl
│
├── scripts/
│   ├── run_1d.jl           # Example: 1D simulation
│   ├── run_2d.jl           # Example: 2D temperature
│   ├── run_2d_moisture.jl  # Example: Full moisture
│   └── run_2d_twolayer.jl  # Example: Circulation
│
└── output/
    └── runs/               # Timestamped results
```

## Physics Model

### Orbital Mechanics

The moon has a 28-hour rotation period (day length) and a 64-hour orbital period around its parent planet. This non-synchronous rotation means:

- Every location experiences day and night as the moon rotates
- The subsolar point moves westward at ~12.9°/hour
- Eclipses occur when the moon passes behind the planet (4-hour duration every 64 hours)

### Energy Balance

The core equation for each grid cell:

```
dT/dt = (Q_in - Q_out + Q_transport + Q_latent) / C
```

Where:
- **Q_in** = Solar input × (1 - albedo) × atmospheric transmission
- **Q_out** = εσT⁴ × IR transmissivity (greenhouse-modified)
- **Q_transport** = Heat diffusion from neighbors (terrain-weighted)
- **Q_latent** = Latent heat from precipitation
- **C** = Biome-dependent heat capacity

### Key Feedbacks

**Ice-Albedo:** Smooth tanh transition at 273K
```julia
ice_fraction = 0.5 * (1 - tanh((T - 273.15) / 10))
albedo = 0.28 + 0.32 * ice_fraction
```

**Water Vapor Greenhouse:** Logarithmic saturation
```julia
τ_water = 0.8 * log(1 + M / 5.0)
transmissivity = 2 / (2 + τ_base + τ_water)
```

**Orographic Precipitation:** Altitude-enhanced condensation
```julia
T_altitude = T - 0.0065 * elevation * 5000  # Lapse rate
M_sat = clausius_clapeyron(T_altitude)
precip = max(0, M - M_sat) * PRECIP_RATE
```

## Configuration

Key parameters in `src/constants.jl`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ROTATION_PERIOD` | 28 hr | Day length |
| `ORBITAL_PERIOD` | 64 hr | Orbit around planet |
| `ECLIPSE_DURATION` | 4 hr | Shadow period per orbit |
| `SOLAR_CONSTANT` | 2000 W/m² | Stellar flux |
| `BASE_HEAT_TRANSPORT` | 15 W/m²/K | Conduction coefficient |
| `CONVECTION_BOOST` | 30 W/m²/K | Added above 270K |
| `BASE_IR_OPTICAL_DEPTH` | 1.0 | Background greenhouse |
| `EVAP_RATE` | 3×10⁻⁶ | Evaporation coefficient |

## Output Files

Simulations save to `output/runs/YYYY-MM-DD_HHMMSS_type/`:

```
run_2025-01-15_143022_2d_moisture/
├── README.md           # Run metadata and statistics
├── solution.jls        # Full ODE solution (JLD2)
├── moon.jls           # Moon configuration
└── plots/
    ├── temperature_final.png
    ├── moisture_final.png
    ├── precipitation.png
    └── biomes.png
```

Load saved results:
```julia
using JLD2
sol = load("output/runs/.../solution.jls", "solution")
moon = load("output/runs/.../moon.jls", "moon")
```

## Performance

| Grid Size | Mode | Time (1000 hr sim) | Threads |
|-----------|------|-------------------|---------|
| 18 lat | 1D | ~2 sec | 1 |
| 45×90 | 2D temp | ~30 sec | 4 |
| 90×180 | 2D + moisture | ~3 min | 4 |
| 45×90 | Two-layer | ~5 min | 4 |

## License

MIT License - See LICENSE file for details.

## Acknowledgments

Built with:
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) - ODE solver
- [Plots.jl](https://github.com/JuliaPlots/Plots.jl) - Visualization
- [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) - Data serialization

## References

- Pierrehumbert, R.T. (2010). *Principles of Planetary Climate*
- Wordsworth, R. (2015). "Atmospheric Heat Redistribution and Collapse on Tidally Locked Rocky Planets"
- Koll, D.D.B. & Abbot, D.S. (2016). "Temperature Structure and Atmospheric Circulation of Dry Tidally Locked Rocky Exoplanets"
