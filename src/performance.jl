"""
Performance utilities: threading control and optimization helpers
"""

"""
    set_threading(enabled::Bool)

Enable or disable multithreading for the simulation.
Threading is only beneficial for larger grids (> THREADING_MIN_CELLS).

When enabled and Julia has multiple threads, hot loops in the solver
will use Threads.@threads for parallelization.
"""
function set_threading(enabled::Bool)
    USE_THREADING[] = enabled && Threads.nthreads() > 1
end

"""
    should_thread(n_cells::Int)

Determine if threading should be used based on grid size.
Small grids have too much overhead to benefit from threading.

Returns true if:
- Threading is enabled (USE_THREADING[] == true)
- Grid has at least THREADING_MIN_CELLS cells
"""
@inline function should_thread(n_cells::Int)
    return USE_THREADING[] && n_cells >= THREADING_MIN_CELLS
end

"""
    get_threading_status()

Return current threading configuration as a named tuple.
Useful for debugging and logging.
"""
function get_threading_status()
    return (
        enabled = USE_THREADING[],
        available_threads = Threads.nthreads(),
        min_cells = THREADING_MIN_CELLS
    )
end
