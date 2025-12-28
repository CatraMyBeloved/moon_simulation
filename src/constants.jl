"""
Physical constants for the Hot Moon simulation
"""

# Stefan-Boltzmann constant (W⋅m⁻²⋅K⁻⁴)
const STEFAN_BOLTZMANN = 5.67e-8

# Solar constant at moon's orbit (W⋅m⁻²)
const SOLAR_CONSTANT = 2200.0

# Orbital parameters (seconds)
const ROTATION_PERIOD = 28.0 * 3600    # 28 hours
const ORBITAL_PERIOD = 64.0 * 3600     # 64 hours
const ECLIPSE_DURATION = 4.0 * 3600    # 4 hours

# Surface properties
const EMISSIVITY = 0.95                # Dimensionless

# Heat transport - non-linear with convection
# At low temps: only conduction (slow). At high temps: convection kicks in (fast).
const HEAT_TRANSPORT_BASE = 10.0       # W⋅m⁻²⋅K⁻¹ (conduction only, always present)
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
const HEAT_CAPACITY_OCEAN = 4.0e7      # Deep water - massive thermal inertia
const HEAT_CAPACITY_WETLAND = 2.0e7    # Coastal/wetland - moist soil
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
const BASE_GREENHOUSE = 0.45
const OPTICAL_DEPTH = 0.25

# Ice-albedo feedback
const ICE_ALBEDO = 0.60                # Increased from 0.50 (fresh snow ~0.8-0.9, sea ice ~0.5-0.7)
const ICE_TRANSITION_WIDTH = 10.0      # K
const FREEZING_POINT = 273.15          # K

# Water vapor feedback
const WATER_VAPOR_REF_TEMP = 273.15     # K
const WATER_VAPOR_SCALE = 30.0         # K
const WATER_VAPOR_STRENGTH = 0.08

# Atmospheric scattering
const SCATTER_FRACTION = 0.35