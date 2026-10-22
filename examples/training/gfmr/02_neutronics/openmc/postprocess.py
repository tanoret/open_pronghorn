"""Post-process the GFMR OpenMC reference: k-eff summary, kinetics, power maps and the 8-group,
6-precursor-group diffusion constants for the MOOSE neutronics (02_neutronics/xs/).

Two homogenizations are written:
  core       volume-weighted fuel salt + graphite of the active lattice (one unit cell is
             ~6.6 % salt by volume); used in the RZ diffusion model for the core region
  reflector  graphite
The core constants are flux-volume weighted from the material-wise MGXS of the nominal case:
  Sigma_x,core,g = (sum_m Sigma_x,m,g phi_m,g V_m) / (sum_m phi_m,g V_m)
which is exact for reaction rates when the intra-cell flux is resolved by the tallies (it is: the
material tallies integrate over the actual heterogeneous geometry). The graphite material MGXS
include both the moderator and the reflector; they are reused for the reflector region.

Usage: python postprocess.py --ref DIR --out ../xs
"""

import argparse
import json
import math
import os
import re
import sys

import numpy as np
import openmc
import openmc.mgxs

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import gfmr  # noqa: E402


def keff(d):
    sp = [f for f in os.listdir(d) if f.startswith('statepoint')]
    with openmc.StatePoint(os.path.join(d, sorted(sp)[-1])) as s:
        return s.keff, s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ref', required=True)
    ap.add_argument('--out', default=os.path.join(HERE, '..', 'xs'))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    summary = {}

    # ------------------------------------------------------------- reactivity coefficients
    ks = {}
    for c in ('nominal', 'rods_out', 'rods_in', 'fuel_hot', 'graphite_hot'):
        d = os.path.join(args.ref, c)
        if os.path.isdir(d) and any(f.startswith('statepoint') for f in os.listdir(d)):
            k, _ = keff(d)
            ks[c] = k
            summary['k_' + c] = [k.nominal_value, k.std_dev]

    def rho_pcm(k):
        return (k - 1.0) / k * 1e5

    if 'rods_out' in ks:
        w = rho_pcm(ks['rods_out']) - rho_pcm(ks['nominal'])
        summary['rod_worth_36cm_pcm'] = [w.nominal_value, w.std_dev]
    if 'rods_in' in ks and 'rods_out' in ks:
        w = rho_pcm(ks['rods_out']) - rho_pcm(ks['rods_in'])
        summary['rod_worth_total_pcm'] = [w.nominal_value, w.std_dev]
    if 'fuel_hot' in ks:
        a = (rho_pcm(ks['fuel_hot']) - rho_pcm(ks['nominal'])) / 100.0
        summary['alpha_fuel_pcm_per_K'] = [a.nominal_value, a.std_dev]
    if 'graphite_hot' in ks:
        a = (rho_pcm(ks['graphite_hot']) - rho_pcm(ks['nominal'])) / 100.0
        summary['alpha_graphite_pcm_per_K'] = [a.nominal_value, a.std_dev]

    # ------------------------------------------------------------- nominal: kinetics, maps
    d = os.path.join(args.ref, 'nominal')
    k, sp = keff(d)
    sp_path = os.path.join(d, sorted(f for f in os.listdir(d) if f.startswith('statepoint'))[-1])
    s = openmc.StatePoint(sp_path)
    kin = s.get_kinetics_parameters()
    summary['generation_time_s'] = [kin.generation_time.nominal_value, kin.generation_time.std_dev]
    summary['beta_eff_total'] = [kin.beta_effective.nominal_value, kin.beta_effective.std_dev]
    tg = s.get_tally(name='ifp_groups')
    tden = s.get_tally(name='ifp')
    den = tden.get_values(scores=['ifp-denominator']).flatten()[0]
    beta_g = tg.get_values(scores=['ifp-beta-numerator']).flatten() / den
    summary['beta_eff_groups'] = beta_g.tolist()

    # energy deposition split. The 'heating' case has coupled neutron-photon transport, so its
    # 'heating' score deposits the prompt / capture gammas where they are absorbed. (In the
    # neutron-only nominal case the photon energy is deposited locally at the collision site.)
    sm = openmc.Summary(os.path.join(d, 'summary.h5'))
    dh = os.path.join(args.ref, 'heating')
    src = dh if os.path.isdir(dh) and any(f.startswith('statepoint') for f in os.listdir(dh)) else d
    sh = openmc.StatePoint(os.path.join(src, sorted(
        f for f in os.listdir(src) if f.startswith('statepoint'))[-1]))
    smh = openmc.Summary(os.path.join(src, 'summary.h5'))
    id2name = {m.id: m.name for m in smh.materials}
    gh = sh.get_tally(name='global_heating')
    heat = gh.get_values(scores=['heating']).flatten()
    summary['heating_source_case'] = os.path.basename(src)
    summary['heating_fraction_by_material'] = {
        id2name[int(b)]: float(h / heat.sum()) for b, h in zip(gh.filters[0].bins, heat)}
    summary['graphite_heating_fraction'] = summary['heating_fraction_by_material']['graphite']

    # per-channel power peaking (kappa-fission, normalized to the mean of the 379 channels).
    # Distribcell instance k of a salt cell has the path '...->l<id>(i,j)->...', where (i, j)
    # are the axial lattice indices of gfmr.lattice().
    sm.geometry.determine_paths()
    all_cells = sm.geometry.get_all_cells()
    rows = []
    for zone, cname in (('inner', 'inner_salt'), ('outer', 'outer_salt'), ('cr', 'cr_salt')):
        tally = s.get_tally(name='channel_power_' + cname)
        cell = next(c for c in all_cells.values() if c.name == cname)
        kf = tally.get_values(scores=['kappa-fission']).flatten()
        for path, v in zip(cell.paths, kf):
            i, j = (int(n) for n in re.search(r'l\d+\((-?\d+),(-?\d+)\)', path).groups())
            x = gfmr.PITCH * (i + j / 2.0)
            y = gfmr.PITCH * math.sqrt(3) / 2.0 * j
            rows.append((zone, i, j, x, y, v))
    kf_mean = np.mean([r[5] for r in rows])
    with open(os.path.join(args.out, 'channel_power.csv'), 'w') as f:
        f.write('type,i,j,x,y,relative_power,channel_power_W\n')
        for t, i, j, x, y, v in rows:
            f.write('%s,%d,%d,%.6f,%.6f,%.6f,%.6e\n' % (
                t, i, j, x, y, v / kf_mean, v / sum(r[5] for r in rows) * gfmr.POWER
                * (1 - summary['graphite_heating_fraction'])))
    summary['channel_peaking'] = float(max(r[5] for r in rows) / kf_mean)
    for zone in ('inner', 'outer', 'cr'):
        summary['zone_power_fraction_' + zone] = float(
            sum(r[5] for r in rows if r[0] == zone) / sum(r[5] for r in rows))

    # axial power shape (from the R-Z flux tally, kappa-fission summed over radius & groups)
    rz = s.get_tally(name='rz_flux')
    # mesh bins: r fastest, then phi, then z; the energy filter varies fastest of all
    kfz = rz.get_values(scores=['kappa-fission']).reshape(36, 1, 35, 8)[:, 0, :, :].sum(axis=2)
    axial = kfz.sum(axis=1)
    axial = axial / axial.mean()
    summary['axial_peaking'] = float(axial.max())
    z = (np.arange(36) + 0.5) * gfmr.HEIGHT / 36
    np.savetxt(os.path.join(args.out, 'axial_power.csv'), np.c_[z, axial], delimiter=',',
               header='z,relative_power', comments='')

    # ------------------------------------------------------------- multigroup constants
    lib = openmc.mgxs.Library.load_from_file('mgxs_library', directory=d)
    lib.load_from_statepoint(s)
    by_name = {c.name: c for c in lib.domains}
    G = lib.energy_groups.num_groups

    def xs(cell, t):
        # fast group first
        return np.asarray(lib.get_mgxs(by_name[cell], t).get_xs(order_groups='increasing'))

    def flux(cell):
        # domain-integrated flux (flux x volume) per group; the energy filter runs low -> high
        m = lib.get_mgxs(by_name[cell], 'total')
        return m.tallies['flux'].get_values().flatten()[::-1]

    def constants(cell):
        """Diffusion constants of one cell region (SI conversion is done downstream)."""
        r = {'flux_volume': flux(cell).tolist()}
        tr = xs(cell, 'transport')
        # void regions (He) have ~zero cross sections: cap D at the region size scale
        r['D'] = (1.0 / (3.0 * np.maximum(tr, 1e-3))).tolist()
        r['sigma_a'] = xs(cell, 'absorption').tolist()
        r['nu_sigma_f'] = xs(cell, 'nu-fission').tolist()
        r['kappa_sigma_f'] = xs(cell, 'kappa-fission').tolist()
        r['sigma_f'] = xs(cell, 'fission').tolist()
        r['sigma_s'] = xs(cell, 'nu-scatter matrix').tolist()      # [g_from][g_to]
        r['inv_v'] = xs(cell, 'inverse-velocity').tolist()
        r['rel_err_absorption_max'] = float(np.nanmax(np.asarray(lib.get_mgxs(
            by_name[cell], 'absorption').get_xs(value='rel_err'))))
        return r

    out = {'groups_eV_decreasing': lib.energy_groups.group_edges[::-1].tolist(),
           'num_groups': G, 'num_precursor_groups': 6, 'regions': {}}
    for name in by_name:
        out['regions'][name] = constants(name)

    # Homogenized lattice for the R-Z model: flux-volume weighting over every region of the
    # active lattice (all cells but the radial reflector); reaction rates are preserved exactly
    lattice = [n for n in by_name if n != 'reflector_lattice']

    def hom(t, matrix=False):
        num, den = 0.0, 0.0
        for n in lattice:
            w = flux(n)
            num = num + xs(n, t) * (w[:, None] if matrix else w)
            den = den + w
        return num / (den[:, None] if matrix else den)

    def hom_D():
        # flux-weighted transport cross section, D = 1 / (3 Sigma_tr)
        num, den = 0.0, 0.0
        for n in lattice:
            w = flux(n)
            num = num + xs(n, 'transport') * w
            den = den + w
        return 1.0 / (3.0 * num / den)

    out['core'] = {'D': hom_D().tolist(), 'sigma_a': hom('absorption').tolist(),
                   'nu_sigma_f': hom('nu-fission').tolist(),
                   'kappa_sigma_f': hom('kappa-fission').tolist(),
                   'sigma_f': hom('fission').tolist(),
                   'sigma_s': hom('nu-scatter matrix', matrix=True).tolist(),
                   'inv_v': hom('inverse-velocity').tolist()}
    out['reflector'] = {k: v for k, v in out['regions']['reflector_lattice'].items()
                        if k in ('D', 'sigma_a', 'nu_sigma_f', 'kappa_sigma_f', 'sigma_f',
                                 'sigma_s', 'inv_v')}

    # spectra and precursor data, from the fuel salt of the inner zone
    out['chi_prompt'] = xs('inner_salt', 'chi-prompt').tolist()
    out['chi_delayed'] = xs('inner_salt', 'chi-delayed').tolist()
    out['beta_material'] = xs('inner_salt', 'beta').tolist()
    out['decay_constants'] = xs('inner_salt', 'decay-rate').tolist()
    out['beta_eff'] = summary['beta_eff_groups']
    out['generation_time'] = summary['generation_time_s'][0]
    with open(os.path.join(args.out, 'gfmr_8g.json'), 'w') as f:
        json.dump(out, f, indent=1)
    with open(os.path.join(args.out, 'openmc_summary.json'), 'w') as f:
        json.dump(summary, f, indent=1)
    for k_, v in summary.items():
        print(k_, v)


if __name__ == '__main__':
    main()
