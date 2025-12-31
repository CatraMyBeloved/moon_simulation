"""
Grid iteration utilities for 2D climate simulations.

Provides abstractions for iterating over grid cells with automatic threading
decisions and neighbor access patterns.
"""

# =============================================================================
# Direction Constants
# =============================================================================

const DIRECTION_NORTH = 1
const DIRECTION_SOUTH = 2
const DIRECTION_EAST = 3
const DIRECTION_WEST = 4

const ALL_DIRECTIONS = (DIRECTION_NORTH, DIRECTION_SOUTH, DIRECTION_EAST, DIRECTION_WEST)

# =============================================================================
# Neighbor Access Functions
# =============================================================================

"""
    get_north_neighbor_indices(i, j, n_lat) -> (exists::Bool, ni::Int, nj::Int)

Return neighbor indices for the cell to the north (lower latitude index).
Returns (false, 0, 0) at the northern boundary.
"""
@inline function get_north_neighbor_indices(i::Int, j::Int, n_lat::Int)
    if i > 1
        return (true, i - 1, j)
    else
        return (false, 0, 0)
    end
end

"""
    get_south_neighbor_indices(i, j, n_lat) -> (exists::Bool, ni::Int, nj::Int)

Return neighbor indices for the cell to the south (higher latitude index).
Returns (false, 0, 0) at the southern boundary.
"""
@inline function get_south_neighbor_indices(i::Int, j::Int, n_lat::Int)
    if i < n_lat
        return (true, i + 1, j)
    else
        return (false, 0, 0)
    end
end

"""
    get_east_neighbor_index(j, n_lon) -> Int

Return longitude index of the eastern neighbor (wraps around).
Always exists due to periodic boundary.
"""
@inline function get_east_neighbor_index(j::Int, n_lon::Int)
    return mod1(j + 1, n_lon)
end

"""
    get_west_neighbor_index(j, n_lon) -> Int

Return longitude index of the western neighbor (wraps around).
Always exists due to periodic boundary.
"""
@inline function get_west_neighbor_index(j::Int, n_lon::Int)
    return mod1(j - 1, n_lon)
end

"""
    get_neighbor_indices(i, j, direction, n_lat, n_lon) -> (exists::Bool, ni::Int, nj::Int)

Return neighbor indices for the specified direction.
For east/west, always returns exists=true due to periodic boundary.
"""
@inline function get_neighbor_indices(i::Int, j::Int, direction::Int, n_lat::Int, n_lon::Int)
    if direction == DIRECTION_NORTH
        return get_north_neighbor_indices(i, j, n_lat)
    elseif direction == DIRECTION_SOUTH
        return get_south_neighbor_indices(i, j, n_lat)
    elseif direction == DIRECTION_EAST
        return (true, i, get_east_neighbor_index(j, n_lon))
    else  # DIRECTION_WEST
        return (true, i, get_west_neighbor_index(j, n_lon))
    end
end

# =============================================================================
# Cell Iteration Functions
# =============================================================================

"""
    foreach_cell(f, n_lat, n_lon)

Execute function f(i, j) for each grid cell, with automatic threading.

Uses threading when the grid has enough cells (>THREADING_MIN_CELLS) and
threading is enabled. Otherwise uses a simple nested loop with SIMD hints.
"""
function foreach_cell(f, n_lat::Int, n_lon::Int)
    n_cells = n_lat * n_lon
    if should_thread(n_cells)
        @inbounds Threads.@threads for idx in 1:n_cells
            i = ((idx - 1) % n_lat) + 1
            j = ((idx - 1) ÷ n_lat) + 1
            f(i, j)
        end
    else
        @inbounds for i in 1:n_lat
            for j in 1:n_lon
                f(i, j)
            end
        end
    end
end

"""
    foreach_cell_indexed(f, n_lat, n_lon)

Execute function f(i, j, idx) for each grid cell, including linear index.

Uses threading when beneficial. The linear index idx ranges from 1 to n_lat*n_lon.
"""
function foreach_cell_indexed(f, n_lat::Int, n_lon::Int)
    n_cells = n_lat * n_lon
    if should_thread(n_cells)
        @inbounds Threads.@threads for idx in 1:n_cells
            i = ((idx - 1) % n_lat) + 1
            j = ((idx - 1) ÷ n_lat) + 1
            f(i, j, idx)
        end
    else
        idx = 0
        @inbounds for i in 1:n_lat
            for j in 1:n_lon
                idx += 1
                f(i, j, idx)
            end
        end
    end
end

"""
    foreach_cell_sequential(f, n_lat, n_lon)

Execute function f(i, j) for each grid cell, without threading.

Use this when the operation has data dependencies that prevent parallelization,
or when threading overhead would outweigh benefits.
"""
function foreach_cell_sequential(f, n_lat::Int, n_lon::Int)
    @inbounds for i in 1:n_lat
        for j in 1:n_lon
            f(i, j)
        end
    end
end

# =============================================================================
# Neighbor Iteration
# =============================================================================

"""
    foreach_valid_neighbor(f, i, j, n_lat, n_lon)

Execute function f(ni, nj, direction) for each valid neighbor of cell (i,j).

Skips neighbors that don't exist (at latitude boundaries).
East/west neighbors always exist due to periodic boundary.
"""
@inline function foreach_valid_neighbor(f, i::Int, j::Int, n_lat::Int, n_lon::Int)
    # North neighbor
    if i > 1
        f(i - 1, j, DIRECTION_NORTH)
    end

    # South neighbor
    if i < n_lat
        f(i + 1, j, DIRECTION_SOUTH)
    end

    # East neighbor (always exists, wraps)
    f(i, mod1(j + 1, n_lon), DIRECTION_EAST)

    # West neighbor (always exists, wraps)
    f(i, mod1(j - 1, n_lon), DIRECTION_WEST)
end

# =============================================================================
# Index Conversion
# =============================================================================

"""
    linear_to_grid_indices(idx, n_lat) -> (i::Int, j::Int)

Convert linear index to (latitude, longitude) grid indices.
"""
@inline function linear_to_grid_indices(idx::Int, n_lat::Int)
    i = ((idx - 1) % n_lat) + 1
    j = ((idx - 1) ÷ n_lat) + 1
    return (i, j)
end

"""
    grid_to_linear_index(i, j, n_lat) -> Int

Convert (latitude, longitude) grid indices to linear index.
"""
@inline function grid_to_linear_index(i::Int, j::Int, n_lat::Int)
    return (j - 1) * n_lat + i
end
