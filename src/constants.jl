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
const EVAP_RATE = 1.0e-5                # Base evaporation rate (kg/m²/s per K above threshold)
const EVAP_THRESHOLD = 280.0            # Minimum temperature for evaporation (K)

# Precipitation
const PRECIP_RATE = 0.001

# Moisture transport
# Note: Unlike heat transport which is divided by heat capacity (~1e6), moisture
# transport directly affects dM/dt. So this value must be much smaller to get
# similar timescales (hours, not milliseconds).
const MOISTURE_DIFFUSION = 1e-4       # Base moisture diffusion coefficient (1/s)
const MOISTURE_BARRIER_STRENGTH = 6.0   # Mountains block moisture more than heat

# Latent heat of vaporization (J/kg)
const LATENT_HEAT = 2.5e6