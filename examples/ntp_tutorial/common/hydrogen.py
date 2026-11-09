"""Ideal-gas hydrogen with NASA polynomials, in Python.

The same equations as the MOOSE fluid properties NTPHydrogenFluidProperties
(``common/hydrogen_nasa.i``), for the post-processing scripts of the tutorial.
``python3 01_closures/check_hydrogen.py`` compares the two.

Model
-----
* thermal ideal gas            rho = p / (R T)
* NASA 7-coefficient fit       cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
  (GRI-Mech 3.0, break at 1000 K, valid up to 3500 K; cp, mu and k are held at their
  T_MAX values above it and h continues linearly, as in NTPHydrogenFluidProperties)
* Sutherland viscosity         mu = 8.76e-6 (T/300)^1.5 (300 + 72) / (T + 72)
* Eucken conductivity          k  = mu (cp + 5/4 R)

There is no dissociation and no ortho/para conversion. The absolute enthalpy
level is that of the NASA fit (h = 0 near 298 K).
"""

from __future__ import annotations

import numpy as np

MOLAR_MASS = 2.01588e-3  # kg/mol
R_UNIVERSAL = 8.314462618  # J/(mol K)
R = R_UNIVERSAL / MOLAR_MASS  # J/(kg K)

_LOW = (2.34433112, 7.98052075e-03, -1.94781510e-05, 2.01572094e-08, -7.37611761e-12, -9.17935173e02)
_HIGH = (3.33727920, -4.94024731e-05, 4.99456778e-07, -1.79566394e-10, 2.00255376e-14, -9.50158922e02)

T_MIN = 20.0
T_MAX = 3500.0


def _clamp(T):
    return np.clip(np.asarray(T, dtype=float), T_MIN, T_MAX)


def _coeffs(T):
    T = np.asarray(T, dtype=float)
    hi = T >= 1000.0
    return T, hi


def cp(T):
    """Isobaric specific heat [J/(kg K)], held constant outside [T_MIN, T_MAX]."""
    T, hi = _coeffs(_clamp(T))
    out = np.empty_like(T)
    for mask, a in ((~hi, _LOW), (hi, _HIGH)):
        t = T[mask]
        out[mask] = a[0] + a[1] * t + a[2] * t**2 + a[3] * t**3 + a[4] * t**4
    return out * R


def h(T):
    """Specific enthalpy [J/kg] on the NASA reference, continued linearly outside [T_MIN, T_MAX]."""
    T_raw = np.asarray(T, dtype=float)
    Tc, hi = _coeffs(_clamp(T_raw))
    out = np.empty_like(Tc)
    for mask, a in ((~hi, _LOW), (hi, _HIGH)):
        t = Tc[mask]
        out[mask] = a[0] * t + a[1] / 2 * t**2 + a[2] / 3 * t**3 + a[3] / 4 * t**4 + a[4] / 5 * t**5 + a[5]
    return out * R + cp(Tc) * (T_raw - Tc)


def e(T):
    """Specific internal energy [J/kg]."""
    return h(T) - R * np.asarray(T, dtype=float)


def mu(T):
    """Dynamic viscosity [Pa s], Sutherland, held constant outside [T_MIN, T_MAX]."""
    T = _clamp(T)
    return 8.76e-6 * (T / 300.0) ** 1.5 * (300.0 + 72.0) / (T + 72.0)


def k(T):
    """Thermal conductivity [W/(m K)], Eucken."""
    return mu(T) * (cp(T) + 1.25 * R)


def rho(T, p):
    """Density [kg/m^3]."""
    return np.asarray(p, dtype=float) / (R * np.asarray(T, dtype=float))


def gamma(T):
    c = cp(T)
    return c / (c - R)


def T_from_h(h_target, T_lo=T_MIN, T_hi=2 * T_MAX):
    """Invert h(T) by bisection; h is monotone (linear beyond T_MAX)."""
    h_target = np.asarray(h_target, dtype=float)
    lo = np.full_like(h_target, T_lo)
    hi = np.full_like(h_target, T_hi)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        below = h(mid) < h_target
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)
    return 0.5 * (lo + hi)
