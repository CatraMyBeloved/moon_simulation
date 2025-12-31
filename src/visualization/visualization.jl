"""
Visualization module: plotting and animation for climate simulations.

This module provides visualization functions for both 1D and 2D simulations.
Key design: generic functions that work with any PlotVariable type to eliminate duplication.
"""

# Include all visualization submodules
include("viz_helpers.jl")
include("viz_1d.jl")
include("viz_2d_snapshots.jl")
include("viz_2d_timeseries.jl")
include("viz_animation.jl")
