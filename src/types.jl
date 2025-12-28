"""
Data structures for Hot Moon simulation

Uses Julia's multiple dispatch pattern: abstract type with concrete 1D/2D implementations.
"""

"""
    AbstractMoonBody

Abstract supertype for all moon body representations.
Enables multiple dispatch for simulation functions.
"""
abstract type AbstractMoonBody end

# =============================================================================
# 1D Model - Latitude bands only
# =============================================================================

"""
    MoonBody1D <: AbstractMoonBody

1D climate model with latitude bands (no longitude variation).
Useful for quick simulations and understanding latitudinal patterns.

# Fields
- `n_lat::Int`: Number of latitude bands (equator to pole)
- `latitudes::Vector{Float64}`: Center latitude of each band (degrees)
- `cell_areas::Vector{Float64}`: Normalized area weights (cos-weighted)
- `_transport_cache::Vector{Float64}`: Pre-allocated buffer for ODE solver
"""
struct MoonBody1D <: AbstractMoonBody
    n_lat::Int
    latitudes::Vector{Float64}
    cell_areas::Vector{Float64}
    _transport_cache::Vector{Float64}
end

"""
    MoonBody1D(n_lat=18)

Create a 1D moon with the specified number of latitude bands.
"""
function MoonBody1D(n_lat::Int=18)
    lat_spacing = 90.0 / n_lat
    latitudes = [(i + 0.5) * lat_spacing for i in 0:n_lat-1]

    # Area weighting: proportional to cos(latitude)
    cell_areas = [cos(deg2rad(lat)) for lat in latitudes]
    cell_areas ./= sum(cell_areas)

    transport_cache = zeros(Float64, n_lat)

    return MoonBody1D(n_lat, latitudes, cell_areas, transport_cache)
end

function Base.show(io::IO, moon::MoonBody1D)
    print(io, "MoonBody1D($(moon.n_lat) latitude bands)")
end

# =============================================================================
# 2D Model - Latitude × Longitude grid
# =============================================================================

"""
    MoonBody2D <: AbstractMoonBody

2D climate model with latitude/longitude grid.
Full spatial representation with topography support.

# Fields
- `n_lat::Int`: Number of latitude bands (pole to pole, -90° to +90°)
- `n_lon::Int`: Number of longitude cells
- `latitudes::Vector{Float64}`: Center latitude of each band (degrees)
- `longitudes::Vector{Float64}`: Center longitude of each cell (degrees)
- `cell_areas::Matrix{Float64}`: Normalized area weights [lat, lon]
- `elevation::Matrix{Float64}`: Surface elevation in meters [lat, lon]
- `transport_coeffs::Array{Float64,3}`: Heat transport modifiers [lat, lon, direction] (1=N, 2=S, 3=E, 4=W)
- `_transport_cache::Vector{Float64}`: Pre-allocated buffer for ODE solver
- `moisture_transport_coeffs::Array{Float64,3}`: Moisture transport modifiers [lat, lon, direction]
- `_moisture_cache::Vector{Float64}`: Pre-allocated buffer for moisture transport
"""
struct MoonBody2D <: AbstractMoonBody
    n_lat::Int
    n_lon::Int
    latitudes::Vector{Float64}
    longitudes::Vector{Float64}
    cell_areas::Matrix{Float64}
    elevation::Matrix{Float64}
    transport_coeffs::Array{Float64,3}
    _transport_cache::Vector{Float64}
    moisture_transport_coeffs::Array{Float64,3}
    _moisture_cache::Vector{Float64}
end

"""
    MoonBody2D(n_lat=18, n_lon=36; seed=42, sea_level=0.1, scale=0.02, octaves=3)

Create a 2D moon with the specified grid resolution.
Generates terrain using fractal noise - elevation < 0 is ocean, > 0 is land.

# Arguments
- `n_lat`: Number of latitude bands (pole to pole, -90° to +90°)
- `n_lon`: Number of longitude cells
- `seed`: Random seed for terrain generation (default: 42)
- `sea_level`: Higher values = more ocean coverage (default: 0.1)
- `scale`: Noise scale - smaller = larger landmasses (default: 0.02)
- `octaves`: Noise detail - fewer = smoother terrain (default: 3)
"""
function MoonBody2D(n_lat::Int=18, n_lon::Int=36; seed::Int=42, sea_level::Float64=0.1,
                    scale::Float64=0.02, octaves::Int=3)
    # Cell centers
    lat_spacing = 90.0 / n_lat
    lon_spacing = 360.0 / n_lon
    latitudes = [(i - 0.5) * 180.0 / n_lat - 90.0 for i in 1:n_lat]
    longitudes = [(j + 0.5) * lon_spacing for j in 0:n_lon-1]

    # Area weighting: proportional to cos(latitude)
    cell_areas = [cos(deg2rad(lat)) for lat in latitudes, _ in 1:n_lon]
    cell_areas ./= sum(cell_areas)

    # Generate terrain using fractal noise
    sampler = create_elevation_sampler(seed=seed, octaves=octaves)
    elevation = [generate_elevation(latitudes[i], longitudes[j], sampler,
                                    sea_level=sea_level, scale=scale)
                 for i in 1:n_lat, j in 1:n_lon]

    # Calculate directional transport coefficients based on elevation
    # Includes: slope barriers, asymmetric downslope flow, ocean bonus
    transport_coeffs = zeros(Float64, n_lat, n_lon, 4)
    calculate_transport_coefficients!(transport_coeffs, elevation, n_lat, n_lon)

    # Calculate moisture transport coefficients (stronger barrier effect)
    moisture_transport_coeffs = zeros(Float64, n_lat, n_lon, 4)
    calculate_moisture_transport_coefficients!(moisture_transport_coeffs, elevation, n_lat, n_lon)

    # Pre-allocated caches for ODE solver
    transport_cache = zeros(Float64, n_lat * n_lon)
    moisture_cache = zeros(Float64, n_lat * n_lon)

    return MoonBody2D(n_lat, n_lon, latitudes, longitudes, cell_areas,
                      elevation, transport_coeffs, transport_cache,
                      moisture_transport_coeffs, moisture_cache)
end

function Base.show(io::IO, moon::MoonBody2D)
    print(io, "MoonBody2D($(moon.n_lat) lat × $(moon.n_lon) lon)")
end

# =============================================================================
# Unified constructor - dispatches based on arguments
# =============================================================================

"""
    HotMoonBody(n_lat) -> MoonBody1D
    HotMoonBody(n_lat, n_lon) -> MoonBody2D

Convenience constructor that returns the appropriate moon type.

# Examples
```julia
moon1d = HotMoonBody(18)        # Returns MoonBody1D
moon2d = HotMoonBody(18, 36)    # Returns MoonBody2D
```
"""
HotMoonBody(n_lat::Int) = MoonBody1D(n_lat)
HotMoonBody(n_lat::Int, n_lon::Int; kwargs...) = MoonBody2D(n_lat, n_lon; kwargs...)
