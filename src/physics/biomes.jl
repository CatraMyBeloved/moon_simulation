"""
Biome classification and heat capacity blending.
"""

# =============================================================================
# Temperature Zone Weights
# =============================================================================

"""
    compute_ice_zone_weight(T) -> Float64

Compute weight for ice sheet zone (T < -10°C).
"""
@inline @fastmath function compute_ice_zone_weight(T::Real)
    return 0.5 * (1 - tanh((T - BIOME_T_ICE) / BIOME_T_WIDTH))
end

"""
    compute_tundra_zone_weight(T) -> Float64

Compute weight for tundra zone (-10°C ≤ T < 0°C).
"""
@inline @fastmath function compute_tundra_zone_weight(T::Real)
    lower = 0.5 * (1 + tanh((T - BIOME_T_ICE) / BIOME_T_WIDTH))
    upper = 0.5 * (1 - tanh((T - BIOME_T_TUNDRA) / BIOME_T_WIDTH))
    return lower * upper
end

"""
    compute_cold_zone_weight(T) -> Float64

Compute weight for cold zone (0°C ≤ T < 10°C).
"""
@inline @fastmath function compute_cold_zone_weight(T::Real)
    lower = 0.5 * (1 + tanh((T - BIOME_T_TUNDRA) / BIOME_T_WIDTH))
    upper = 0.5 * (1 - tanh((T - BIOME_T_COLD) / BIOME_T_WIDTH))
    return lower * upper
end

"""
    compute_temperate_zone_weight(T) -> Float64

Compute weight for temperate zone (10°C ≤ T < 20°C).
"""
@inline @fastmath function compute_temperate_zone_weight(T::Real)
    lower = 0.5 * (1 + tanh((T - BIOME_T_COLD) / BIOME_T_WIDTH))
    upper = 0.5 * (1 - tanh((T - BIOME_T_TEMPERATE) / BIOME_T_WIDTH))
    return lower * upper
end

"""
    compute_hot_zone_weight(T) -> Float64

Compute weight for hot zone (T ≥ 20°C).
"""
@inline @fastmath function compute_hot_zone_weight(T::Real)
    return 0.5 * (1 + tanh((T - BIOME_T_TEMPERATE) / BIOME_T_WIDTH))
end

"""
    compute_all_temperature_zone_weights(T) -> NamedTuple

Compute all temperature zone weights at once.
"""
@inline function compute_all_temperature_zone_weights(T::Real)
    return (
        ice = compute_ice_zone_weight(T),
        tundra = compute_tundra_zone_weight(T),
        cold = compute_cold_zone_weight(T),
        temperate = compute_temperate_zone_weight(T),
        hot = compute_hot_zone_weight(T)
    )
end

# =============================================================================
# Moisture Weights
# =============================================================================

"""
    compute_dry_moisture_weight(M) -> Float64

Compute weight for dry conditions (M < BIOME_M_DRY).
"""
@inline @fastmath function compute_dry_moisture_weight(M::Real)
    return 0.5 * (1 - tanh((M - BIOME_M_DRY) / BIOME_M_WIDTH))
end

"""
    compute_wet_moisture_weight(M) -> Float64

Compute weight for wet conditions (M ≥ BIOME_M_MODERATE).
"""
@inline @fastmath function compute_wet_moisture_weight(M::Real)
    return 0.5 * (1 + tanh((M - BIOME_M_MODERATE) / BIOME_M_WIDTH))
end

"""
    compute_desert_moisture_weight(M) -> Float64

Compute weight for desert conditions (very dry).
"""
@inline @fastmath function compute_desert_moisture_weight(M::Real)
    return 0.5 * (1 - tanh((M - BIOME_M_DESERT) / BIOME_M_WIDTH))
end

"""
    compute_forest_moisture_weight(M) -> Float64

Compute weight for forest conditions (very wet).
"""
@inline @fastmath function compute_forest_moisture_weight(M::Real)
    return 0.5 * (1 + tanh((M - BIOME_M_WET) / BIOME_M_WIDTH))
end

# =============================================================================
# Elevation Weights
# =============================================================================

"""
    compute_mountain_elevation_weight(elevation) -> Float64

Compute weight for mountain terrain (elevation > BIOME_ELEV_MOUNTAIN).
"""
@inline @fastmath function compute_mountain_elevation_weight(elevation::Real)
    return 0.5 * (1 + tanh((elevation - BIOME_ELEV_MOUNTAIN) / BIOME_ELEV_WIDTH))
end

"""
    compute_wetland_elevation_weight(elevation) -> Float64

Compute weight for wetland terrain (0 < elevation < BIOME_ELEV_WETLAND).
"""
@inline @fastmath function compute_wetland_elevation_weight(elevation::Real)
    return 0.5 * (1 - tanh((elevation - BIOME_ELEV_WETLAND) / BIOME_ELEV_WIDTH))
end

# =============================================================================
# Heat Capacity Blending
# =============================================================================

"""
    blend_cold_zone_heat_capacity(wet_weight) -> Float64

Blend heat capacity for cold zone between boreal forest and cold steppe.
"""
@inline function blend_cold_zone_heat_capacity(wet_weight::Real)
    return wet_weight * BIOME_HEAT_CAPACITY[BIOME_BOREAL_FOREST + 1] +
           (1 - wet_weight) * BIOME_HEAT_CAPACITY[BIOME_COLD_STEPPE + 1]
end

"""
    blend_temperate_zone_heat_capacity(wet_weight) -> Float64

Blend heat capacity for temperate zone between forest and grassland.
"""
@inline function blend_temperate_zone_heat_capacity(wet_weight::Real)
    return wet_weight * BIOME_HEAT_CAPACITY[BIOME_TEMPERATE_FOREST + 1] +
           (1 - wet_weight) * BIOME_HEAT_CAPACITY[BIOME_GRASSLAND + 1]
end

"""
    blend_hot_zone_heat_capacity(M) -> Float64

Blend heat capacity for hot zone between tropical forest, savanna, and desert.
"""
@inline function blend_hot_zone_heat_capacity(M::Real)
    desert_weight = compute_desert_moisture_weight(M)
    forest_weight = compute_forest_moisture_weight(M)
    savanna_weight = max(0.0, 1.0 - desert_weight - forest_weight)

    return forest_weight * BIOME_HEAT_CAPACITY[BIOME_TROPICAL_FOREST + 1] +
           savanna_weight * BIOME_HEAT_CAPACITY[BIOME_SAVANNA + 1] +
           desert_weight * BIOME_HEAT_CAPACITY[BIOME_HOT_DESERT + 1]
end

"""
    blend_climate_zone_heat_capacity(temp_weights, M) -> Float64

Blend heat capacity across all temperature zones for a given moisture level.
"""
@inline function blend_climate_zone_heat_capacity(temp_weights, M::Real)
    wet_weight = compute_wet_moisture_weight(M)

    hc_cold = blend_cold_zone_heat_capacity(wet_weight)
    hc_temperate = blend_temperate_zone_heat_capacity(wet_weight)
    hc_hot = blend_hot_zone_heat_capacity(M)

    return temp_weights.ice * BIOME_HEAT_CAPACITY[BIOME_ICE_SHEET + 1] +
           temp_weights.tundra * BIOME_HEAT_CAPACITY[BIOME_TUNDRA + 1] +
           temp_weights.cold * hc_cold +
           temp_weights.temperate * hc_temperate +
           temp_weights.hot * hc_hot
end

"""
    compute_biome_blended_heat_capacity(T, M, elevation) -> Float64

Calculate heat capacity using smooth biome blending.
Handles ocean, mountain, wetland, and climate-based biomes with tanh transitions.
"""
@inline @fastmath function compute_biome_blended_heat_capacity(T::Real, M::Real, elevation::Real)
    # Ocean has depth-scaled heat capacity
    if elevation < 0.0
        depth_factor = clamp(1.0 - elevation, 1.0, 2.0)
        return BIOME_HEAT_CAPACITY[BIOME_OCEAN + 1] * depth_factor
    end

    # Mountain check (early exit for high elevations)
    mountain_weight = compute_mountain_elevation_weight(elevation)
    if mountain_weight > 0.99
        return BIOME_HEAT_CAPACITY[BIOME_MOUNTAIN + 1]
    end

    # Wetland weight for coastal areas
    wetland_weight = compute_wetland_elevation_weight(elevation)

    # Climate-based heat capacity
    temp_weights = compute_all_temperature_zone_weights(T)
    hc_climate = blend_climate_zone_heat_capacity(temp_weights, M)

    # Apply wetland blending
    hc_land = wetland_weight * BIOME_HEAT_CAPACITY[BIOME_WETLAND + 1] +
              (1 - wetland_weight) * hc_climate

    # Apply mountain blending
    return mountain_weight * BIOME_HEAT_CAPACITY[BIOME_MOUNTAIN + 1] +
           (1 - mountain_weight) * hc_land
end

# =============================================================================
# Biome Classification (for visualization)
# =============================================================================

"""
    classify_cell_biome(T, M, elevation) -> Int

Classify a cell into one of 12 biome types for visualization.
Returns biome ID (0-11).
"""
function classify_cell_biome(T::Real, M::Real, elevation::Real)
    # Elevation-based biomes take priority
    if elevation < 0.0
        return BIOME_OCEAN
    elseif elevation > BIOME_ELEV_MOUNTAIN
        return BIOME_MOUNTAIN
    elseif elevation < BIOME_ELEV_WETLAND
        return BIOME_WETLAND
    end

    # Temperature-based classification
    if T < BIOME_T_ICE
        return BIOME_ICE_SHEET
    elseif T < BIOME_T_TUNDRA
        return BIOME_TUNDRA
    elseif T < BIOME_T_COLD
        return M >= BIOME_M_DRY ? BIOME_BOREAL_FOREST : BIOME_COLD_STEPPE
    elseif T < BIOME_T_TEMPERATE
        return M >= BIOME_M_MODERATE ? BIOME_TEMPERATE_FOREST : BIOME_GRASSLAND
    else
        if M >= BIOME_M_WET
            return BIOME_TROPICAL_FOREST
        elseif M >= BIOME_M_DESERT
            return BIOME_SAVANNA
        else
            return BIOME_HOT_DESERT
        end
    end
end

# =============================================================================
# Legacy Heat Capacity Functions (latitude-based, for 1D and simple cases)
# =============================================================================

"""
    compute_latitude_based_heat_capacity(lat_deg) -> Float64

Return heat capacity based on latitude zone.
Used for 1D simulations and as fallback.
"""
@inline function compute_latitude_based_heat_capacity(lat_deg::Real)
    abs_lat = abs(lat_deg)
    if abs_lat < ZONE_EQUATORIAL_END
        return HEAT_CAPACITY_EQUATORIAL
    elseif abs_lat < ZONE_MIDLAT_END
        return HEAT_CAPACITY_MIDLAT
    else
        return HEAT_CAPACITY_POLAR
    end
end

"""
    compute_elevation_based_heat_capacity(lat_deg, elevation) -> Float64

Return heat capacity based on elevation and latitude.
Used for 2D temperature-only simulations.
"""
@inline function compute_elevation_based_heat_capacity(lat_deg::Real, elevation::Real)
    # Ocean
    if elevation < 0.0
        depth_factor = clamp(1.0 - elevation, 1.0, 2.0)
        return HEAT_CAPACITY_OCEAN * depth_factor
    end

    land_cap = compute_latitude_based_heat_capacity(lat_deg)

    # Wetlands zone
    if elevation < ELEV_WETLAND_END
        t = elevation / ELEV_WETLAND_END
        return HEAT_CAPACITY_WETLAND * (1 - t) + land_cap * t
    end

    # Mountain zone
    if elevation > ELEV_MOUNTAIN_START
        t = clamp((elevation - ELEV_MOUNTAIN_START) / (1.0 - ELEV_MOUNTAIN_START), 0.0, 1.0)
        return land_cap * (1 - t) + HEAT_CAPACITY_MOUNTAIN * t
    end

    return land_cap
end
