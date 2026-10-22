"""GFMR reference specification and derived quantities shared by every stage of the training.

Source: Saltfoss Energy ApS, "Generic FUNaK-Fuelled Thermal Molten Salt Reactor (GFMR) Model",
https://saltfossenergy.github.io/GFMR/ (Apache 2.0, DOI 10.5281/zenodo.23348200). Values marked
SPEC are copied from the specification; values marked ASSUMED are not given there and are chosen
here (each is discussed in the README); DERIVED values follow from the two.

The channel map is not tabulated in the specification. It was decoded from the published top
view (geometry-xy.png) by sampling the figure colour at every lattice position, and it reproduces
the stated totals exactly: 82 inner + 288 outer + 9 control-rod channels = 379 fuel channels.
"""

import math

# ----------------------------------------------------------------------------------------------
# SPEC: main design parameters
# ----------------------------------------------------------------------------------------------
POWER = 150.0e6            # W thermal
HEIGHT = 3.60              # m, core (graphite) height
T_FUEL_MEAN_C = 650.0      # degC, fuel temperature (taken as the core-average salt temperature)
T_GRAPHITE_C = 800.0       # degC, graphite temperature (neutronics reference)
MDOT = 2300.0              # kg/s total
ENRICHMENT = 2.12          # wt% U-235 (OpenMC enrichment parameter)
R_REFLECTOR = 1.75         # m, reflector outer radius (3.5 m "core diameter")
D_CHANNEL = 0.040          # m, fuel channel diameter
PITCH = 0.140              # m, hexagon pitch (flat-to-flat of the graphite unit cell)
CR_INSERTION = 0.36        # m, nominal control-rod insertion from the top

# control-rod channel
R_ABSORBER = 0.0065        # m, Hf absorber radius (1.3 cm diameter)
T_GAP = 0.001              # m, helium gap
T_GUIDE = 0.002            # m, Hastelloy N guide tube
R_GAP_OUT = R_ABSORBER + T_GAP            # 0.0075 m
R_GUIDE_OUT = R_GAP_OUT + T_GUIDE         # 0.0095 m

# hydraulics (SPEC)
U_INNER = 1.70             # m/s
U_OUTER = 1.03             # m/s
U_CR = 1.74                # m/s, in the annulus between guide tube and channel wall
DP_INNER = 14.3e3          # Pa
DP_OUTER = 6.0e3           # Pa
DP_CR = 34.0e3             # Pa

# ----------------------------------------------------------------------------------------------
# SPEC: materials
# ----------------------------------------------------------------------------------------------
def rho_salt(T_C):
    """FUNaK density, kg/m^3 (SPEC: 4.808 - 0.00113 T g/cm^3, T in degC)."""
    return (4.808 - 0.00113 * T_C) * 1000.0


def mu_salt(T_C):
    """FUNaK dynamic viscosity, Pa s (SPEC: 736.58 exp(-0.006 T) mPa s, T in degC)."""
    return 736.58e-3 * math.exp(-0.006 * T_C)


CP_SALT = 900.0            # J/kg/K (SPEC)
# ASSUMED: the specification gives no conductivity ("awaiting specification"). 0.8 W/m/K is
# representative of NaF-KF-UF4 / NaF-ZrF4-UF4 fuel salts (ORNL-4541 class data, 0.6-1.0 W/m/K).
K_SALT = 0.8               # W/m/K
K_GRAPHITE = 25.0          # W/m/K (SPEC, constant)
RHO_GRAPHITE = 1700.0      # kg/m^3 (SPEC)
# ASSUMED: graphite heat capacity at ~800 C (nuclear graphite, ~1.9-2.0 kJ/kg/K)
CP_GRAPHITE = 1900.0       # J/kg/K

# salt composition (SPEC): 26.4 UF4 - 24.7 KF - 48.9 NaF mol%
X_UF4, X_KF, X_NAF = 0.264, 0.247, 0.489
M_U = 238.03 - ENRICHMENT / 100 * 3.0     # g/mol, approximate (U-235/U-238 mix)
M_SALT = X_UF4 * (M_U + 4 * 18.998) + X_KF * (39.098 + 18.998) + X_NAF * (22.990 + 18.998)

# ----------------------------------------------------------------------------------------------
# Lattice (DERIVED from the SPEC figure). Axial hex coordinates (i, j): x = p (i + j/2),
# y = p sqrt(3)/2 j, i.e. flat-topped rows along x as in the published top view.
# ring(i, j) = max(|i|, |j|, |i + j|).
#   rings 0-5  : inner zone (91 positions, of which 9 hold control rods -> 82 inner channels)
#   rings 6-10 : outer zone, complete
#   ring 11    : outer zone without the 3 positions nearest each of the 6 corners (48 of 66)
#   control rods: (1,-1), (-1,0), (0,1) [ring 1] and the 6 corners of ring 4.
# ----------------------------------------------------------------------------------------------
N_RINGS = 12               # rings 0..11


def ring(i, j):
    return max(abs(i), abs(j), abs(i + j))


CR_POSITIONS = [(1, -1), (-1, 0), (0, 1), (4, 0), (-4, 0), (0, 4), (0, -4), (4, -4), (-4, 4)]


def _corner_distance(i, j):
    """Hex-grid distance from (i, j) on ring 11 to the nearest ring corner."""
    corners = [(11, 0), (-11, 0), (0, 11), (0, -11), (11, -11), (-11, 11)]
    return min(max(abs(i - a), abs(j - b), abs((i - a) + (j - b))) for a, b in corners)


def channel_type(i, j):
    """'inner', 'outer', 'cr', or None (graphite / reflector, no channel)."""
    r = ring(i, j)
    if r > 11:
        return None
    if r == 11 and _corner_distance(i, j) <= 1:
        # the corner and the two neighbours on the ring are cut out (3 per corner)
        return None
    if (i, j) in CR_POSITIONS:
        return 'cr'
    return 'inner' if r <= 5 else 'outer'


def lattice():
    """List of (i, j, x, y, type) for every channel position."""
    out = []
    for j in range(-11, 12):
        for i in range(-11, 12):
            t = channel_type(i, j)
            if t is not None:
                out.append((i, j, PITCH * (i + j / 2.0), PITCH * math.sqrt(3) / 2.0 * j, t))
    return out


def counts():
    c = {'inner': 0, 'outer': 0, 'cr': 0}
    for *_, t in lattice():
        c[t] += 1
    return c


# ----------------------------------------------------------------------------------------------
# DERIVED hydraulics
# ----------------------------------------------------------------------------------------------
A_CHANNEL = math.pi / 4 * D_CHANNEL**2
A_CR_ANNULUS = math.pi * ((D_CHANNEL / 2) ** 2 - R_GUIDE_OUT**2)
DH_CR = D_CHANNEL - 2 * R_GUIDE_OUT                  # annulus hydraulic diameter
A_CELL = math.sqrt(3) / 2 * PITCH**2                 # graphite unit-cell area


def mass_flow_from_velocities(T_C=T_FUEL_MEAN_C):
    c = counts()
    rho = rho_salt(T_C)
    return rho * (c['inner'] * U_INNER * A_CHANNEL + c['outer'] * U_OUTER * A_CHANNEL
                  + c['cr'] * U_CR * A_CR_ANNULUS)


def delta_T_core(T_C=T_FUEL_MEAN_C):
    return POWER / (MDOT * CP_SALT)


DT_CORE = POWER / (MDOT * CP_SALT)                   # 72.5 K
T_IN_C = T_FUEL_MEAN_C - DT_CORE / 2                 # ASSUMED: 650 C is the core-average salt T
T_OUT_C = T_FUEL_MEAN_C + DT_CORE / 2


def blasius_dp(u, D, L=HEIGHT, T_C=T_FUEL_MEAN_C):
    """Smooth-pipe Darcy friction loss over L (Blasius / McAdams, turbulent)."""
    rho, mu = rho_salt(T_C), mu_salt(T_C)
    Re = rho * u * D / mu
    f = 0.316 * Re**-0.25 if Re < 2e4 else 0.184 * Re**-0.2
    return f * L / D * 0.5 * rho * u**2, Re


if __name__ == '__main__':
    c = counts()
    print('channels', c, 'total', sum(c.values()))
    print('rho %.1f kg/m3, mu %.4f Pa s, Pr %.1f at %.0f C' % (
        rho_salt(T_FUEL_MEAN_C), mu_salt(T_FUEL_MEAN_C),
        CP_SALT * mu_salt(T_FUEL_MEAN_C) / K_SALT, T_FUEL_MEAN_C))
    print('mdot from the zone velocities %.0f kg/s (SPEC 2300)' % mass_flow_from_velocities())
    print('core dT %.1f K, T_in %.1f C, T_out %.1f C' % (DT_CORE, T_IN_C, T_OUT_C))
    for name, u, D, dp in (('inner', U_INNER, D_CHANNEL, DP_INNER),
                           ('outer', U_OUTER, D_CHANNEL, DP_OUTER),
                           ('cr', U_CR, DH_CR, DP_CR)):
        dpf, Re = blasius_dp(u, D)
        print('%-5s Re %7.0f  dP_friction %5.1f kPa  (SPEC %4.1f kPa)' % (name, Re, dpf / 1e3, dp / 1e3))
