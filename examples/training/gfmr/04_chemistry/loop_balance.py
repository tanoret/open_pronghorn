"""0D loop balances of the GFMR fuel-salt chemistry: inlet values for the 3-D core model and the
verification targets of the 3-D result.

The GFMR specification defines only the core. The circulating salt also spends time outside the
core (pipes, pump, heat exchanger). The ex-core inventory is ASSUMED here: the core salt volume
(1.71 m^3, 3.0 s core transit) plus an ex-core volume of 2x the core salt, typical of graphite
MSR designs (MSRE: core ~ 1/3 of the 2 m^3 fuel salt). The resulting loop time is 9.0 s.

Each tracked quantity c obeys, averaged over the loop,
    dc/dt = S_core * (V_core / V_loop) - (lambda + r_removal) c      [loop average]
and in the once-through core
    c_out - c_in = S_core * t_core     (no removal inside the 3 s core transit)
so the inlet value is the loop-average value minus half of the core increment.

Usage: python loop_balance.py   (prints a table and writes inlet_values.i for the MOOSE input)
"""

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
import gfmr  # noqa: E402

N_A = 6.02214076e23
EV = 1.602176634e-19

# ---------------------------------------------------------------- loop (ASSUMED ex-core part)
rho = gfmr.rho_salt(gfmr.T_FUEL_MEAN_C)
c = gfmr.counts()
V_core = gfmr.HEIGHT * ((c['inner'] + c['outer']) * gfmr.A_CHANNEL + c['cr'] * gfmr.A_CR_ANNULUS)
t_core = V_core * rho / gfmr.MDOT
EXCORE_RATIO = 2.0
V_loop = V_core * (1 + EXCORE_RATIO)
t_loop = V_loop * rho / gfmr.MDOT
f_core = V_core / V_loop

# ---------------------------------------------------------------- power and fission rate
P_salt = gfmr.POWER * (1 - 0.0787)       # W deposited in the salt (OpenMC coupled n-gamma)
E_FISSION = 200.0e6 * EV                 # J per fission, as in [FissionProductInventory]
F_rate = gfmr.POWER / E_FISSION          # fissions/s (all recoverable energy per fission)

# ---------------------------------------------------------------- off-gas (ASSUMED)
# Helium sparging in the pump bowl. The specification has none; MSRE stripped ~ 10-15 % of the
# dissolved xenon per pass in its spray ring. Assumed: 10 % of the loop flow's Xe per loop pass.
STRIP_PER_PASS = 0.10
r_strip = -math.log(1 - STRIP_PER_PASS) / t_loop       # 1/s, loop-average removal rate

# ---------------------------------------------------------------- thermal flux for Xe-135 burnout
# one-group thermal flux in the salt, from the OpenMC-derived diffusion of stage 2: the salt
# absorbs mostly thermal neutrons; phi_th = fission rate / (Sigma_f,th V_salt) with the
# thermal-group fission cross section of the salt
xs = json.load(open(os.path.join(HERE, '..', '02_neutronics', 'xs', 'gfmr_8g.json')))
salt = xs['regions']['inner_salt']
flux = salt['flux_volume']                               # fast group first
sig_f = [s * 100.0 for s in salt['sigma_f']]             # 1/m
fis_per_flux = sum(sf * f for sf, f in zip(sig_f, flux)) / sum(flux)
phi_core = F_rate / (fis_per_flux * V_core)              # group-averaged flux in the salt, 1/m^2/s
th = [4, 5, 6, 7]                                        # groups below 4 eV
phi_th_core = phi_core * sum(flux[g] for g in th) / sum(flux)
SIG_XE = 2.6476e6 * 1e-28                                # m^2, database thermal capture (2.65 Mb)

# ---------------------------------------------------------------- redox (U(IV)/U(III) buffer)
# FUNaK: 26.4 mol% UF4. Molar concentration of uranium in the salt:
c_salt = rho / (gfmr.M_SALT * 1e-3)                       # mol salt / m^3
U_total = gfmr.X_UF4 * c_salt
F_total = (4 * gfmr.X_UF4 + gfmr.X_KF + gfmr.X_NAF) * c_salt
X_UIII = 0.01                                             # U(III)/U target (MSRE practice ~1 %, ASSUMED)
EQ_PER_FISSION = 0.8                                      # fission oxidation (ORNL MSRE experience)
ox_rate = EQ_PER_FISSION * F_rate / N_A                   # mol eq / s, whole loop
# redox control: a reductant addition that holds U(III)/U at the target, loop-time constant
TAU_REDOX = 86400.0                                       # s, 1 day (ASSUMED)


def nuclide(name):
    d = json.load(open(os.path.join(HERE, '..', '..', '..', '..', 'data',
                                    'fission_product_database.json')))['nuclides'][name]
    lam = math.log(2) / d['half_life']
    y = d['yields']['U235_thermal']['cumulative']
    return lam, y


def main():
    lam_I, y_I = nuclide('I135')
    lam_Xe, y_Xe = nuclide('Xe135')
    # loop-average concentrations (atoms/m^3), production confined to the core volume
    burn = SIG_XE * phi_th_core * f_core                  # loop-average burnout rate
    I_avg = y_I * F_rate / V_loop / lam_I
    # Xe: direct yield beyond the I-135 cumulative share plus I decay; removal by decay, burnout,
    # stripping
    y_Xe_direct = y_Xe - y_I
    Xe_avg = (y_Xe_direct * F_rate / V_loop + lam_I * I_avg) / (lam_Xe + burn + r_strip)
    Xe_core_rise = y_Xe_direct * F_rate / V_core * t_core
    Xe_in = Xe_avg - 0.5 * Xe_core_rise * f_core
    I_in = I_avg - 0.5 * y_I * F_rate / V_core * t_core * f_core

    # reactivity worth of the xenon (rough, one-group): Sigma_a,Xe / Sigma_a,salt in the thermal range
    sa_th = sum(salt['sigma_a'][g] * 100.0 * flux[g] for g in th) / sum(flux[g] for g in th)
    xe_poison = Xe_avg * SIG_XE / sa_th

    U_III_target = X_UIII * U_total
    # inlet U(III): the core oxidizes during the 3 s transit; the inlet is the loop average plus
    # half of the core oxidation
    dUIII_core = ox_rate / V_core * t_core                # mol/m^3 oxidized per core pass
    U_III_in = U_III_target + 0.5 * dUIII_core

    rows = [('core salt volume', V_core, 'm^3'), ('core transit', t_core, 's'),
            ('loop time (ASSUMED ex-core 2x core)', t_loop, 's'),
            ('fission rate', F_rate, '1/s'), ('thermal flux in the salt', phi_th_core, '1/m^2/s'),
            ('U total', U_total, 'mol/m^3'), ('F total', F_total, 'mol/m^3'),
            ('U(III) target (1 %)', U_III_target, 'mol/m^3'),
            ('fission oxidation, whole loop', ox_rate * 86400, 'mol eq/day'),
            ('U(III) oxidized per core pass', dUIII_core, 'mol/m^3'),
            ('U(III) at the core inlet', U_III_in, 'mol/m^3'),
            ('I-135 loop average', I_avg, '1/m^3'), ('Xe-135 loop average', Xe_avg, '1/m^3'),
            ('Xe-135 stripping rate', r_strip, '1/s'),
            ('Xe-135 burnout rate (loop avg)', burn, '1/s'),
            ('Xe-135 poisoning, Sig_a,Xe / Sig_a,salt (thermal)', xe_poison, '-')]
    for name, v, u in rows:
        print('%-50s %12.4g %s' % (name, v, u))

    with open(os.path.join(HERE, 'inlet_values.i'), 'w') as f:
        f.write('# Generated by loop_balance.py -- do not edit.\n')
        f.write('t_core = %.6g\nt_loop = %.6g\n' % (t_core, t_loop))
        f.write('U_total = %.8g\nF_total = %.8g\n' % (U_total, F_total))
        f.write('U_III_in = %.8g\nU_IV_in = %.8g\n' % (U_III_in, U_total - U_III_in))
        # [FissionProductInventory] works in mol/m^3: convert from atoms/m^3
        f.write('I135_in = %.8g\nXe135_in = %.8g\n' % (I_in / N_A, Xe_in / N_A))
        f.write('I135_loop = %.8g\nXe135_loop = %.8g\n' % (I_avg / N_A, Xe_avg / N_A))
        f.write('phi_thermal = %.8g\n' % phi_th_core)
        f.write('fission_eq_per_J = %.8g\n' % (EQ_PER_FISSION / (E_FISSION * N_A)))
    print('wrote inlet_values.i')


if __name__ == '__main__':
    main()
