"""
Physical constants for the Hot Moon simulation
"""

# Stefan-Boltzmann constant (W⋅m⁻²⋅K⁻⁴)
const STEFAN_BOLTZMANN = 5.67e-8

# Solar constant at moon's orbit (W⋅m⁻²)
const SOLAR_CONSTANT = 2000.0

# Orbital parameters (seconds)
const ROTATION_PERIOD = 28.0 * 3600    # 28 hours
const ORBITAL_PERIOD = 64.0 * 3600     # 64 hours
const ECLIPSE_DURATION = 4.0 * 3600    # 4 hours

# Surface properties
const EMISSIVITY = 0.95                # Dimensionless

# Heat transport - non-linear with convection
# At low temps: only conduction (slow). At high temps: convection kicks in (fast).
const HEAT_TRANSPORT_BASE = 15.0       # W⋅m⁻²⋅K⁻¹ (conduction only, always present)
const HEAT_TRANSPORT_CONVECTION = 30.0 # W⋅m⁻²⋅K⁻¹ (additional convection, activates when warm)
const CONVECTION_THRESHOLD = 270.0     # K (temperature where convection ramps up)
const CONVECTION_WIDTH = 20.0          # K (transition width for tanh)

# Zone-based heat capacity (J⋅m⁻²⋅K⁻¹)
# Heat capacity = ρ × c × d where ρ=density, c=specific heat, d=depth
# These values assume different surface compositions by latitude zone:
const HEAT_CAPACITY_EQUATORIAL = 1.0e6
const HEAT_CAPACITY_MIDLAT = 4.0e6
const HEAT_CAPACITY_POLAR = 1.0e7

# Elevation-based heat capacity
const HEAT_CAPACITY_OCEAN = 2.0e7      # Deep water - massive thermal inertia
const HEAT_CAPACITY_WETLAND = 1.0e7    # Coastal/wetland - moist soil
const HEAT_CAPACITY_MOUNTAIN = 3.0e5   # Bare rock - responds quickly

# Elevation thresholds for terrain types (noise units, roughly -1 to 1 after sea level)
const ELEV_WETLAND_END = 0.15          # Wetlands transition to normal terrain
const ELEV_MOUNTAIN_START = 0.5        # Mountains begin here

# Elevation-based transport modifiers
const SLOPE_BARRIER_STRENGTH = 4.0     # How much steep slopes block transport (higher = more blocking)
const DOWNSLOPE_STRENGTH = 0.4         # Asymmetry: how much easier downhill flow is (0 = symmetric)
const SLOPE_SCALE = 3.0                # Sensitivity of asymmetry to elevation difference
const OCEAN_TRANSPORT_BONUS = 1.5      # Ocean cells transport heat more efficiently (currents)

# Zone boundaries (degrees latitude)
const ZONE_EQUATORIAL_END = 30.0
const ZONE_MIDLAT_END = 70.0

# Atmospheric parameters
const BASE_ALBEDO = 0.28
const SOLAR_OPTICAL_DEPTH = 0.25       # Visible/UV optical depth for solar transmission

# Ice-albedo feedback
const ICE_ALBEDO = 0.60                # Increased from 0.50 (fresh snow ~0.8-0.9, sea ice ~0.5-0.7)
const ICE_TRANSITION_WIDTH = 10.0      # K
const FREEZING_POINT = 273.15          # K

# Infrared optical depth (greenhouse effect)
# Uses Eddington approximation: transmissivity = 2/(2+τ)
# τ=1.0 gives ~33% greenhouse effect, τ=2.0 gives 50%
const BASE_IR_OPTICAL_DEPTH = 1.0      # From permanent gases (CO₂-equivalent)
const WV_OPTICAL_SCALE = 0.8           # Water vapor contribution strength
const WV_REF_MOISTURE = 5.0            # Reference moisture for log scaling (kg/m²)

# Elevation-based cooling (thinner atmosphere at altitude)
# Greenhouse effect is reduced by this fraction at max elevation (elev=1.0)
# This gives roughly 3-4°C/km cooling (real lapse rate is ~6.5°C/km)
const ELEVATION_GREENHOUSE_REDUCTION = 0.5

# ============================================================================
# MOISTURE SYSTEM
# ============================================================================

# Clausius-Clapeyron saturation
const MOISTURE_REF_TEMP = 273.15        # Reference temperature (K)
const MOISTURE_REF_SATURATION = 5.0     # Reference saturation (kg/m²) at 273K
const CLAUSIUS_CLAPEYRON_SCALE = 17.0   # Exponential scale factor

# Lapse rate for orographic lifting
const LAPSE_RATE = 6.5e-3               # Temperature drop per meter (K/m)
const ELEVATION_SCALE = 5000.0          # Max elevation in meters (for normalizing)

# Evaporation
# Earth average: ~4 mm/day ≈ 5×10⁻⁵ kg/m²/s at typical ocean temps
# This moon receives ~1.6× Earth's solar flux, so ~3-4× evaporation is plausible
# At 300K: evap = 1×10⁻⁵ × 20 = 2×10⁻⁴ kg/m²/s ≈ 17 mm/day (4× Earth)
const EVAP_RATE = 0.3e-5                # Base evaporation rate (kg/m²/s per K above threshold)
const EVAP_THRESHOLD = 290.0            # Minimum temperature for evaporation (K)

# Precipitation
const PRECIP_RATE = 0.002
const UPPER_LAYER_PRECIP_RATE = 0.005    # Upper layer precipitation rate (faster than surface - cold air can't hold moisture)
const UPPER_LAYER_TEMP_DROP = 40.0       # K - temperature drop from surface to upper layer for saturation calc

# Moisture transport
# Note: Unlike heat transport which is divided by heat capacity (~1e6), moisture
# transport directly affects dM/dt. So this value must be much smaller to get
# similar timescales (hours, not milliseconds).
const MOISTURE_DIFFUSION = 1e-4         # Base moisture diffusion coefficient (1/s)
const MOISTURE_BARRIER_STRENGTH = 8.0   # Mountains block moisture more than heat
const MOISTURE_DOWNSLOPE_PREFERENCE = 0.3  # Asymmetry factor for downhill moisture flow
const MOISTURE_SLOPE_SENSITIVITY = 3.0  # Sensitivity of moisture flow to slope

# Latent heat of vaporization (J/kg)
const LATENT_HEAT = 2.5e6

# ============================================================================
# BIOME SYSTEM
# ============================================================================

# Biome IDs (0-11)
const BIOME_OCEAN = 0
const BIOME_ICE_SHEET = 1
const BIOME_TUNDRA = 2
const BIOME_BOREAL_FOREST = 3
const BIOME_COLD_STEPPE = 4
const BIOME_TEMPERATE_FOREST = 5
const BIOME_GRASSLAND = 6
const BIOME_TROPICAL_FOREST = 7
const BIOME_SAVANNA = 8
const BIOME_HOT_DESERT = 9
const BIOME_MOUNTAIN = 10
const BIOME_WETLAND = 11

const NUM_BIOMES = 12

# Biome heat capacities (J⋅m⁻²⋅K⁻¹)
# Index is biome_id + 1 (Julia is 1-indexed)
const BIOME_HEAT_CAPACITY = [
    2.0e7,   # 0: Ocean - water's high specific heat
    1.5e7,   # 1: Ice Sheet - thick ice mass
    1.0e7,   # 2: Tundra - permafrost
    8.0e6,   # 3: Boreal Forest - conifers + wet soil
    3.0e6,   # 4: Cold Steppe - sparse grass, dry soil
    6.0e6,   # 5: Temperate Forest - deciduous + moist soil
    4.0e6,   # 6: Grassland - grass + soil
    8.0e6,   # 7: Tropical Forest - dense canopy, very wet
    3.0e6,   # 8: Savanna - sparse trees, dry soil
    5.0e5,   # 9: Hot Desert - sand/bare rock
    3.0e5,   # 10: Mountain - bare rock
    1.2e7,   # 11: Wetland - coastal marsh, saturated soil
]

# Biome names for visualization
const BIOME_NAMES = [
    "Ocean", "Ice Sheet", "Tundra", "Boreal Forest", "Cold Steppe",
    "Temperate Forest", "Grassland", "Tropical Forest", "Savanna",
    "Hot Desert", "Mountain", "Wetland"
]

# Temperature thresholds (Kelvin)
const BIOME_T_ICE = 263.15       # -10°C - below this is ice sheet
const BIOME_T_TUNDRA = 273.15    # 0°C - below this is tundra
const BIOME_T_COLD = 283.15      # 10°C - below this is cold biomes
const BIOME_T_TEMPERATE = 293.15 # 20°C - below this is temperate biomes

# Moisture thresholds (kg/m²) for biome classification
const BIOME_M_DESERT = 2.0       # Below this: desert
const BIOME_M_DRY = 3.0          # Below this: steppe/dry
const BIOME_M_MODERATE = 4.0     # Below this: grassland
const BIOME_M_WET = 5.0          # Above this: forest

# Elevation thresholds for terrain-based biomes
const BIOME_ELEV_MOUNTAIN = 0.5  # Above this: mountain
const BIOME_ELEV_WETLAND = 0.1   # Below this (but above 0): wetland

# Transition widths for smooth blending (tanh)
const BIOME_T_WIDTH = 5.0        # K - temperature transition width
const BIOME_M_WIDTH = 1.0        # kg/m² - moisture transition width
const BIOME_ELEV_WIDTH = 0.05    # elevation transition width

# Initial state estimation parameters
const INIT_T_EQUATOR = 310.0     # K (~37°C) - equilibrium temp at equator
const INIT_T_POLE = 250.0        # K (~-23°C) - equilibrium temp at poles
const INIT_MOISTURE_DIFFUSE_ITERS = 15  # iterations for moisture diffusion
const INIT_MOISTURE_DECAY = 0.7  # decay factor per diffusion step
const INIT_MOISTURE_SATURATION_CAP = 0.8  # cap initial moisture at this fraction of saturation
const DEFAULT_MOISTURE_ESTIMATE_KG_M2 = 3.0  # moisture estimate for temperature-only sims

# ============================================================================
# PERFORMANCE TUNING
# ============================================================================

# Minimum grid cells before enabling threading (overhead not worth it for small grids)
const THREADING_MIN_CELLS = 100

# Enable threading by default when Julia has multiple threads available
const USE_THREADING = Ref(Threads.nthreads() > 1)

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