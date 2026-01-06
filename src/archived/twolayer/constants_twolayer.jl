"""
Two-Layer Atmosphere Constants (ARCHIVED)

This file contains constants for the two-layer atmosphere implementation
that has been archived due to parameter tuning challenges.

To restore: include this file and the corresponding solver/atmosphere modules.
"""

# Upper layer precipitation (was in main constants)
const UPPER_LAYER_PRECIP_RATE = 0.005    # Upper layer precipitation rate (faster than surface - cold air can't hold moisture)
const UPPER_LAYER_TEMP_DROP = 40.0       # K - temperature drop from surface to upper layer for saturation calc

# ============================================================================
# TWO-LAYER ATMOSPHERE
# ============================================================================

# Ascent parameters
const THERMAL_RESPONSE_SCALE = 20.0      # K - sensitivity to temperature anomaly
const MIN_MOISTURE_CONVECTION = 2.0      # kg/m² - minimum M for deep convection (lowered for easier triggering)
const ASCENT_RATE_MAX = 5e-5             # 1/s - maximum ascent rate (halved to reduce moisture overshoot)

# Absolute convection trigger (replaces relative zonal mean comparison)
const CONVECTION_TEMP_THRESHOLD = 295.0  # K (~22°C) - minimum surface T for convection
const MOISTURE_ASCENT_SENSITIVITY = 0.5  # log-scale multiplier for moisture boost

# Descent parameters
const BASE_DESCENT_RATE = 5e-6           # 1/s - background subsidence everywhere (halved for stability)
const MASS_DESCENT_COEFF = 1.5e-4        # 1/s per unit excess U (reduced slightly)
const POLAR_SINK_STRENGTH = 2e-5         # 1/s - enhanced descent at poles (reduced)
const POLAR_SINK_LAT = 60.0              # degrees - where polar sink begins

# Moisture during ascent
const LIFT_TEMPERATURE_DROP = 30.0       # K - adiabatic cooling during lift
const MOISTURE_SURVIVE_FRACTION = 0.3    # fraction surviving to upper layer
const LIFT_FRACTION = 0.1                # fraction of surface moisture lifted per unit ascent

# Upper layer transport
const UPPER_MERIDIONAL_COEFF = 2.5e-4    # 1/s - poleward mass flow rate (raised from 1e-4)
const UPPER_ZONAL_COEFF = 5e-5           # 1/s - zonal pressure flow
const WESTERLY_BIAS_STRENGTH = 2e-5      # 1/s - Coriolis-like westerly tendency

# Descent suppression of precipitation
const DESCENT_DRYING_SCALE = 2.0         # multiplicative factor on saturation per unit descent

# Numerical safety
const U_FLOOR = 0.3                      # minimum U for divisions (raised from 0.1)
const U_FLOOR_RESTORATION_RATE = 0.5     # rate of restoration when U < U_FLOOR (raised from 0.1)
const U_INITIAL = 1.0                    # initial upper mass (uniform)
const M_UP_INITIAL = 0.05                # kg/m² - initial upper moisture

# ============================================================================
# UPPER LAYER TEMPERATURE (Phase 2 - Full Two-Layer Atmosphere)
# ============================================================================

# Initial conditions
const T_UP_INITIAL = 250.0              # K - initial upper temperature (~-23°C)
const T_UP_EQUILIBRIUM = 240.0          # K - radiative equilibrium (warmer due to 2000 W/m² solar)

# Heat capacity and radiation
const UPPER_LAYER_HEAT_CAPACITY = 1e7   # J/m²/K - effective heat capacity of upper layer
const UPPER_RADIATIVE_COOLING_RATE = 2e-6  # 1/s - radiative relaxation rate

# Vertical heat exchange
const ASCENT_HEAT_TRANSFER_FRACTION = 0.3  # fraction of ΔT transferred during ascent
const DESCENT_ADIABATIC_WARMING = 30.0     # K - warming during full descent

# Temperature-circulation feedback
const VERTICAL_INSTABILITY_SCALE = 50.0    # K - sensitivity of ascent to T_surf - T_up
const T_DESCENT_SENSITIVITY = 0.5          # descent enhancement per K cold anomaly

# Numerical stability for T_up
const T_UP_FLOOR = 180.0                # K - minimum upper temperature
const T_UP_CEILING = 320.0              # K - maximum upper temperature
