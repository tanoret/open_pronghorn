"""Nusselt number of salt_pipe.i against the Gnielinski correlation.

Reads salt_pipe_out.csv (wall temperature, wall heat flux) and the radial profile at z = 1.61 m
(L/D = 40) and compares Nu = q'' D / (k (T_w - T_b)) with
  Nu_G = (f/8)(Re - 1000) Pr / (1 + 12.7 sqrt(f/8)(Pr^(2/3) - 1)),  f = (0.790 ln Re - 1.64)^-2
(Gnielinski 1976, valid 0.5 < Pr < 2000, 3e3 < Re < 5e6, quoted uncertainty +-10-20 %).

Usage: python nusselt.py [run directory]
"""

import csv
import glob
import math
import os
import sys

d = sys.argv[1] if len(sys.argv) > 1 else '.'
rho, mu, cp, k = 4073.5, 0.01494, 900.0, 0.8
D, R, R_out = 0.04, 0.02, 0.035

pp = list(csv.DictReader(open(os.path.join(d, 'salt_pipe_out.csv'))))[-1]
prof = sorted(glob.glob(os.path.join(d, 'salt_pipe_out_radial_0*.csv')))[-1]
rows = list(csv.DictReader(open(prof)))
r = [float(x['x']) for x in rows]
T = [float(x['T_fluid']) for x in rows]
v = [float(x['vel_y']) for x in rows]
edges = [0.0] + [(r[i] + r[i + 1]) / 2 for i in range(len(r) - 1)] + [R]
w = [edges[i + 1] ** 2 - edges[i] ** 2 for i in range(len(r))]
Tb = sum(T[i] * v[i] * w[i] for i in range(len(r))) / sum(v[i] * w[i] for i in range(len(r)))
ub = sum(v[i] * w[i] for i in range(len(r))) / sum(w)

q_wall = float(pp['heat_in']) / (2 * math.pi * R * 2.0)  # uniform heating -> uniform flux
# wall temperature: from the first solid cell with the exact conduction profile of a uniformly
# heated annulus insulated at R_out (the graphite sleeve), T(r) - T(R) =
#   q_s/(4 k_s) (R^2 - r^2) + q_s R_out^2/(2 k_s) ln(r/R)
k_s = 25.0
q_s = float(pp['heat_in']) / (math.pi * (R_out**2 - R**2) * 2.0)
sol = sorted(glob.glob(os.path.join(d, 'salt_pipe_out_radial_solid_*.csv')))[-1]
srow = list(csv.DictReader(open(sol)))[0]
r1, T1 = float(srow['x']), float(srow['T_solid'])
Tw = T1 - (q_s / (4 * k_s) * (R**2 - r1**2) + q_s * R_out**2 / (2 * k_s) * math.log(r1 / R))
Nu = q_wall * D / (k * (Tw - Tb))

Re = rho * ub * D / mu
Pr = cp * mu / k
f = (0.790 * math.log(Re) - 1.64) ** -2
Nu_G = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * math.sqrt(f / 8) * (Pr ** (2 / 3) - 1))
# Dittus-Boelter (heating) for reference
Nu_DB = 0.023 * Re**0.8 * Pr**0.4

print('Re %.0f  Pr %.2f  u_b %.3f m/s' % (Re, Pr, ub))
print('T_wall %.3f K  T_bulk %.3f K  q %.0f W/m2' % (Tw, Tb, q_wall))
print('Nu CFD %.1f   Gnielinski %.1f (ratio %.3f)   Dittus-Boelter %.1f' % (
    Nu, Nu_G, Nu / Nu_G, Nu_DB))
y1 = R - r[-1]
f_p = (0.790 * math.log(Re) - 1.64) ** -2
print('first cell y %.2e m  y+ (Petukhov u_tau) %.1f' % (y1, y1 * ub * math.sqrt(f_p / 8) * rho / mu))
print('energy balance %.2e' % float(pp['energy_balance']))
