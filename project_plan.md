## Project: Hot Moon Climate Simulator

### Overview
Build a 2D climate simulation for a moon with mountains, using energy balance ODEs (no fluid dynamics). Visualize temperature patterns, day/night cycles, and topographic effects.

---

## Project Structure

```
hotmoon/
├── src/
│   ├── types.jl              # Data structures (HotMoon, HotMoon2D)
│   ├── constants.jl          # Physical constants
│   ├── physics.jl            # Physics functions (albedo, greenhouse, solar)
│   ├── geometry.jl           # Sun position, zenith angles
│   ├── atmosphere.jl         # Atmospheric effects (transmission, scattering)
│   ├── topography.jl         # Mountain generation, transport coefficients
│   ├── solver_1d.jl          # 1D simulation (latitude bands only)
│   ├── solver_2d.jl          # 2D simulation (lat + lon grid)
│   └── visualization.jl      # Plotting and animation
│
├── notebooks/
│   ├── exploration.jl        # Pluto notebook for parameter exploration
│   └── analysis.jl           # Pluto notebook for analyzing results
│
├── scripts/
│   ├── run_1d.jl            # Run 1D simulation
│   ├── run_2d.jl            # Run 2D simulation
│   └── compare_scenarios.jl  # Compare different configurations
│
├── test/
│   ├── test_physics.jl      # Unit tests for physics functions
│   └── test_geometry.jl     # Unit tests for geometry
│
├── output/
│   ├── plots/               # Generated plots
│   ├── animations/          # GIFs and videos
│   └── data/                # Saved simulation results
│
├── docs/
│   ├── physics_notes.md     # Physics equations and references
│   └── parameters.md        # Parameter choices and justification
│
├── Project.toml             # Julia dependencies
└── README.md                # Project overview
```

---

## Development Phases

### **Phase 0: Setup (30 min)**
**Goal:** Get environment ready

```julia
# Create project
mkdir hotmoon && cd hotmoon
julia --project=.

# Add packages
using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Plots")
Pkg.add("PlutoUI")
Pkg.add("Pluto")
```

**Files to create:**
- `Project.toml` (auto-generated)
- `README.md` (basic description)

---

### **Phase 1: 1D Foundation (4-6 hours)**
**Goal:** Working 1D simulation with placeholder physics

#### Step 1.1: Core structure (1 hour)
**Create:** `src/constants.jl`
```julia
const STEFAN_BOLTZMANN = 5.67e-8
const SOLAR_CONSTANT = 1361.0
const ROTATION_PERIOD = 28.0 * 3600
const ORBITAL_PERIOD = 64.0 * 3600
const ECLIPSE_DURATION = 4.0 * 3600
const HEAT_CAPACITY = 2.0e6
const HEAT_TRANSPORT = 40.0
const EMISSIVITY = 0.90
```

**Create:** `src/types.jl`
```julia
struct HotMoon
    n_bands::Int
    latitudes::Vector{Float64}
    band_areas::Vector{Float64}
end

function HotMoon(n_bands::Int=18)
    # Constructor implementation
end
```

**Create:** `src/physics.jl`
```julia
# Placeholder implementations
get_albedo(moon::HotMoon, T::Float64) = 0.3
get_greenhouse(moon::HotMoon, T::Float64) = 0.6
# ... etc
```

**Create:** `src/solver_1d.jl`
```julia
function derivatives_1d!(dT, temps, moon::HotMoon, t)
    # ODE system
end

function run_simulation_1d(moon, hours, T0)
    # Solve and return results
end
```

**Create:** `scripts/run_1d.jl`
```julia
include("../src/constants.jl")
include("../src/types.jl")
include("../src/physics.jl")
include("../src/solver_1d.jl")

moon = HotMoon(18)
sol = run_simulation_1d(moon, 200.0, fill(285.0, 18))
# Basic plot
```

**Milestone:** Run `julia scripts/run_1d.jl` → see plot (even with placeholders)

#### Step 1.2: Real physics (2-3 hours)
**Implement in `src/physics.jl`:**
- [ ] `get_albedo` - ice-albedo feedback (tanh function)
- [ ] `get_greenhouse` - water vapor feedback (tanh function)
- [ ] Look up equations online as you go

**Implement in `src/geometry.jl`:**
- [ ] `get_zenith` - sun position calculation
- [ ] `is_eclipsed` - eclipse timing

**Implement in `src/atmosphere.jl`:**
- [ ] `get_solar` - solar heating with atmospheric transmission

**Create:** `test/test_physics.jl`
```julia
using Test
include("../src/physics.jl")

@testset "Physics" begin
    @test get_albedo(moon, 250.0) > get_albedo(moon, 300.0)
    @test get_greenhouse(moon, 290.0) > get_greenhouse(moon, 250.0)
end
```

**Milestone:** Realistic temperature evolution, stable equilibrium

#### Step 1.3: Visualization (1 hour)
**Create:** `src/visualization.jl`
```julia
function plot_global_mean(sol, moon)
    # Time series of global mean
end

function plot_latitude_profile(sol, moon)
    # Mean temperature by latitude
end

function plot_heatmap(sol, moon)
    # Temperature evolution over time
end
```

**Milestone:** Publication-quality plots

---

### **Phase 2: 2D Expansion (6-8 hours)**
**Goal:** Add longitude dimension, proper day/night

#### Step 2.1: 2D structure (2 hours)
**Create:** `src/types.jl` (add to existing)
```julia
struct HotMoon2D
    n_lat::Int
    n_lon::Int
    latitudes::Vector{Float64}
    longitudes::Vector{Float64}
    band_areas::Vector{Float64}
end

function HotMoon2D(n_lat=18, n_lon=36)
    # Constructor
end
```

**Create:** `src/solver_2d.jl`
```julia
function derivatives_2d!(dT, temps, moon::HotMoon2D, t)
    # Reshape to 2D grid
    # Calculate heating/cooling for each grid cell
    # Handle wrapping in longitude
end

function run_simulation_2d(moon, hours, T0)
    # Solve and return
end
```

**Milestone:** 2D grid working (even with simple physics)

#### Step 2.2: Proper solar geometry (2 hours)
**Update:** `src/geometry.jl`
```julia
function get_zenith_2d(moon, t, i, j)
    # Calculate zenith angle for grid cell (i,j)
    # Account for subsolar point movement
end

function get_solar_2d(moon, t, i, j, T)
    # Solar heating for specific grid cell
end
```

**Milestone:** See day/night terminator move across planet

#### Step 2.3: 2D visualization (2 hours)
**Update:** `src/visualization.jl`
```julia
function plot_2d_snapshot(temps, moon)
    # Heatmap of temperature
end

function animate_2d(sol, moon, output_path)
    # Create rotating globe animation
end
```

**Create:** `scripts/run_2d.jl`

**Milestone:** Animated GIF of temperature patterns

---

### **Phase 3: Mountains (3-4 hours)**
**Goal:** Add topography and variable heat transport

#### Step 3.1: Topography generation (1 hour)
**Create:** `src/topography.jl`
```julia
function generate_mountains(n_lat, n_lon; mountain_ranges=1)
    # Create elevation matrix
    # Add mountain ranges at specified locations
end

function create_island_chain(elevation, center_lat, center_lon)
    # Add volcanic islands
end
```

**Update:** `src/types.jl`
```julia
struct HotMoon2D
    # ... existing fields
    elevation::Matrix{Float64}
    transport_coeffs::Array{Float64, 3}  # [lat, lon, direction]
end
```

**Milestone:** Elevation map created

#### Step 3.2: Transport modification (2 hours)
**Add to:** `src/topography.jl`
```julia
function calculate_transport_coefficients(elevation, base_k=40.0)
    # Reduce transport over mountains
    # Direction-dependent (harder to cross mountains)
end
```

**Update:** `src/solver_2d.jl`
- Use `moon.transport_coeffs` instead of constant transport

**Milestone:** See temperature discontinuities at mountain ranges

#### Step 3.3: Analysis (1 hour)
**Create:** `scripts/compare_scenarios.jl`
```julia
# Compare flat vs mountainous
sol_flat = run_simulation_2d(HotMoon2D(18, 36, flat=true), ...)
sol_mountains = run_simulation_2d(HotMoon2D(18, 36, flat=false), ...)

plot_comparison(sol_flat, sol_mountains)
```

**Milestone:** Show mountain climate effects

---

### **Phase 4: Polish & Analysis (2-4 hours)**
**Goal:** Make it pretty and scientifically interesting

#### Tasks:
- [ ] Add climate zone classification
- [ ] Calculate statistics (day-night contrast, pole-equator gradient)
- [ ] Create summary plots
- [ ] Write physics documentation
- [ ] Add parameter sensitivity analysis
- [ ] Create demo Pluto notebook

**Create:** `notebooks/exploration.jl` (Pluto)
```julia
@bind n_bands Slider(10:30, default=18)
@bind rotation_hours Slider(20:40, default=28)
@bind mountain_height Slider(0:500:5000, default=3000)

# Run simulation with these parameters
# Show results interactively
```

---

## Timeline Estimate

### **Tomorrow (8-hour drive as passenger):**
- Hours 1-2: Setup + Phase 1.1 (structure with placeholders)
- Hours 3-5: Phase 1.2 (implement real physics - look up equations)
- Hours 6-7: Phase 1.3 (visualization) + debugging
- Hour 8: Play in Pluto, explore parameters

**End of day:** Working 1D simulation with realistic physics

### **Week 2 (evenings/weekend):**
- Session 1: Phase 2.1-2.2 (2D grid + solar geometry)
- Session 2: Phase 2.3 (2D visualization)
- Session 3: Phase 3.1-3.2 (mountains)
- Session 4: Phase 3.3 + Phase 4 (analysis & polish)

**End of week:** Complete 2D simulation with mountains

### **Week 3+ (optional):**
- Orbital mechanics integration
- More sophisticated features (ice caps, oceans, etc.)
- Publication-quality writeup

---

## Quick Start Commands

```bash
# Day 1 morning:
mkdir hotmoon && cd hotmoon
julia --project=. -e 'using Pkg; Pkg.add(["DifferentialEquations", "Plots", "Pluto", "PlutoUI"])'
mkdir -p src scripts notebooks output/{plots,animations,data} test docs

# Start coding!
julia --project=. scripts/run_1d.jl

# Or start Pluto:
julia --project=. -e 'using Pluto; Pluto.run()'
```

---

## Success Criteria

**Phase 1 Complete When:**
- ✓ 1D simulation runs without errors
- ✓ Temperature reaches stable equilibrium
- ✓ Global mean is reasonable (~0-30°C)
- ✓ Poles are colder than equator
- ✓ Plot looks smooth and physical

**Phase 2 Complete When:**
- ✓ 2D simulation runs without errors
- ✓ Clear day/night terminator visible
- ✓ Hot spot follows sun
- ✓ Animation shows rotating temperature pattern
- ✓ Day-night contrast is reasonable (10-30°C)

**Phase 3 Complete When:**
- ✓ Mountains create temperature discontinuities
- ✓ "Rain shadow" effect visible
- ✓ Comparison plot shows mountain impact
- ✓ Transport coefficients make physical sense

**Phase 4 Complete When:**
- ✓ Can explain every physical effect in the model
- ✓ Parameters are justified
- ✓ Multiple scenarios explored
- ✓ Results are beautiful and interesting
- ✓ Friends say "wow that's cool!"

---

## Key Design Decisions

**Resolution:**
- Start: 18 lat × 36 lon = 648 ODEs (~10° resolution)
- If too slow: reduce to 12 lat × 24 lon = 288 ODEs
- If fast enough: increase to 36 lat × 72 lon = 2,592 ODEs

**Simulation time:**
- Development: 100-200 hours (fast iteration)
- Final runs: 1,000-8,000 hours (full equilibration)

**Output:**
- Save every 1 hour of simulation time
- Downsample for animations (every 10 hours)

**Modularity:**
- Each physics function is independent
- Can swap implementations easily
- Easy to test each piece

