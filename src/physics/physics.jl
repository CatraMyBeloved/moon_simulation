"""
Physics module: collects all physics submodules.

This file includes all physics components and exports the public API.
The original monolithic physics.jl has been split into focused modules:
- radiative.jl: Albedo, greenhouse effect, IR transmissivity
- terrain.jl: Elevation generation, transport coefficients
- moisture.jl: Saturation, evaporation, precipitation
- biomes.jl: Biome classification, heat capacity blending
- atmosphere.jl: Two-layer vertical exchange
- initial_state.jl: Initial condition estimation
"""

# Include all physics submodules
include("radiative.jl")
include("terrain.jl")
include("moisture.jl")
include("biomes.jl")
include("atmosphere.jl")
include("initial_state.jl")
