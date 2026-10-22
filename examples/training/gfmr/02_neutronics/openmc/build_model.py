"""OpenMC model of the GFMR (Saltfoss reference specification), full 3-D core.

Materials are the OpenMC definitions published in the specification (FUNaK, graphite with 1 appm
B, Hf, Hastelloy N, He). Geometry: 379 hexagonal unit cells of 14 cm pitch, 4 cm channels, nine
control-rod channels with the rod inserted 36 cm from the top, a cylindrical graphite reflector to
R = 1.75 m, and the core height of 3.6 m. The specification gives no axial reflector; the top and
bottom planes are vacuum boundaries ("bare" axially), as in the reference description.

Temperatures (SPEC): fuel 650 C, graphite 800 C.

Usage (from this directory, with the GFMR tools environment):
  python build_model.py [--rods in|out] [--particles N] [--batches B] [--inactive I]
                        [--fuel-T C] [--graphite-T C] [--tallies] [--mgxs] [--out DIR]
"""

import argparse
import math
import os
import sys

import numpy as np
import openmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import gfmr  # noqa: E402


def materials(T_fuel_C, T_gr_C, enrichment=gfmr.ENRICHMENT):
    T_fuel_K = T_fuel_C + 273.15
    T_gr_K = T_gr_C + 273.15

    fuel = openmc.Material(name='fuel')
    fuel.add_element('Na', 0.1751432664756447, 'ao')
    fuel.add_element('K', 0.08846704871060172, 'ao')
    fuel.add_element('F', 0.6418338108882521, 'ao')
    fuel.add_element('U', 0.094555874, 'ao', enrichment=enrichment)
    fuel.set_density('g/cm3', 4.808 - 0.00113 * T_fuel_C)
    fuel.temperature = T_fuel_K
    fuel.depletable = True

    moderator = openmc.Material(name='graphite')
    moderator.add_element('C', 0.999999)
    moderator.add_element('B', 0.000001)
    moderator.set_density('g/cm3', 1.7)
    moderator.add_s_alpha_beta('c_Graphite')
    moderator.temperature = T_gr_K

    absorber = openmc.Material(name='control_rod')
    absorber.add_element('Hf', 1.0, 'ao')
    absorber.set_density('g/cm3', 13.3)
    absorber.temperature = T_fuel_K

    guide = openmc.Material(name='hastelloy_n')
    for el, w in (('Ni', 71), ('Cr', 7), ('Mo', 16), ('Fe', 4), ('Si', 1), ('Mn', 0.8), ('V', 0.2)):
        guide.add_element(el, w, percent_type='wo')
    guide.set_density('g/cm3', 8.86)
    guide.temperature = T_fuel_K

    helium = openmc.Material(name='helium')
    helium.add_element('He', 1.0)
    helium.set_density('g/cm3', 0.0001785)
    helium.temperature = T_fuel_K

    return dict(fuel=fuel, graphite=moderator, absorber=absorber, guide=guide, helium=helium)


def geometry(mats, rods_in=True, insertion=None):
    """GFMR core as a hexagonal lattice of five unit-cell universes (inner / outer channel,
    control-rod channel, graphite plug, outside). Cells of 14 cm pitch with flats facing the
    six neighbours (OpenMC 'x' lattice orientation: rows along x, as in the published top view
    and in gfmr.lattice()). Returns (geometry, dict of tally objects)."""
    H = gfmr.HEIGHT * 100.0
    # rod insertion depth from the top (m); default: the nominal 36 cm of the specification
    ins = gfmr.CR_INSERTION if insertion is None else insertion
    z_tip = H - (ins * 100.0 if rods_in else 0.0)
    p = gfmr.PITCH * 100.0
    r_ch = gfmr.D_CHANNEL / 2 * 100.0
    r_abs, r_gap, r_gd = gfmr.R_ABSORBER * 100, gfmr.R_GAP_OUT * 100, gfmr.R_GUIDE_OUT * 100

    cyl = openmc.ZCylinder(r=r_ch)
    c_abs = openmc.ZCylinder(r=r_abs)
    c_gap = openmc.ZCylinder(r=r_gap)
    c_gd = openmc.ZCylinder(r=r_gd)
    z_rod = openmc.ZPlane(z0=z_tip)

    def channel_universe(name):
        salt = openmc.Cell(name=name + '_salt', fill=mats['fuel'], region=-cyl)
        gr = openmc.Cell(name=name + '_graphite', fill=mats['graphite'], region=+cyl)
        return openmc.Universe(name=name, cells=[salt, gr]), salt, gr

    u_inner, salt_inner, gr_inner = channel_universe('inner')
    u_outer, salt_outer, gr_outer = channel_universe('outer')

    salt_cr = openmc.Cell(name='cr_salt', fill=mats['fuel'], region=+c_gd & -cyl)
    gr_cr = openmc.Cell(name='cr_graphite', fill=mats['graphite'], region=+cyl)
    tube = openmc.Cell(name='cr_guide', fill=mats['guide'], region=+c_gap & -c_gd)
    cr_cells = [salt_cr, gr_cr, tube]
    if rods_in:
        cr_cells += [openmc.Cell(name='cr_rod', fill=mats['absorber'], region=-c_abs & +z_rod),
                     openmc.Cell(name='cr_he_low', fill=mats['helium'], region=-c_abs & -z_rod),
                     openmc.Cell(name='cr_gap', fill=mats['helium'], region=+c_abs & -c_gap)]
    else:
        cr_cells.append(openmc.Cell(name='cr_he', fill=mats['helium'], region=-c_gap))
    u_cr = openmc.Universe(name='cr', cells=cr_cells)

    u_plug = openmc.Universe(name='plug', cells=[
        openmc.Cell(name='plug_graphite', fill=mats['graphite'])])
    # outside the 397 lattice positions (still inside the reflector cylinder): graphite
    u_outside = openmc.Universe(name='outside', cells=[
        openmc.Cell(name='reflector_lattice', fill=mats['graphite'])])

    lat = openmc.HexLattice(name='core')
    lat.orientation = 'x'
    lat.center = (0.0, 0.0)
    lat.pitch = (p,)
    n_rings = 13   # one ring of reflector graphite beyond ring 11
    lat.outer = u_outside
    # gfmr uses axial coordinates x = p (i + j/2), y = p sqrt(3)/2 j, which is exactly the
    # (x, alpha) skewed basis of an 'x'-orientation HexLattice (get_local_coordinates), so the
    # cell (i, j) sits at HexLattice.get_universe_index((i, j)).
    univ = {'inner': u_inner, 'outer': u_outer, 'cr': u_cr, None: u_plug}
    rings = [[u_outside] * max(6 * (n_rings - 1 - r), 1) for r in range(n_rings)]
    lat._num_rings = n_rings
    for j in range(-11, 12):
        for i in range(-11, 12):
            if gfmr.ring(i, j) <= 11:
                ring, pos = lat.get_universe_index((i, j))
                rings[ring][pos] = univ[gfmr.channel_type(i, j)]
    lat.universes = rings

    zb = openmc.ZPlane(z0=0.0, boundary_type='vacuum')
    zt = openmc.ZPlane(z0=H, boundary_type='vacuum')
    r_refl = openmc.ZCylinder(r=gfmr.R_REFLECTOR * 100.0, boundary_type='vacuum')
    core = openmc.Cell(name='core', fill=lat, region=-r_refl & +zb & -zt)
    geom = openmc.Geometry([core])
    return geom, dict(salt=[salt_inner, salt_outer, salt_cr], lattice=lat,
                      graphite=[gr_inner, gr_outer, gr_cr])


def settings(particles, batches, inactive):
    s = openmc.Settings()
    s.batches = batches
    s.inactive = inactive
    s.particles = particles
    H = gfmr.HEIGHT * 100.0
    R = 1.40 * 100.0
    s.source = openmc.IndependentSource(
        space=openmc.stats.CylindricalIndependent(
            r=openmc.stats.PowerLaw(0.0, R, 1.0), phi=openmc.stats.Uniform(0, 2 * math.pi),
            z=openmc.stats.Uniform(0.0, H)),
        constraints={'fissionable': True})
    # the salt is ~6.6 % of the lattice volume, so only ~7 % of the uniformly sampled initial
    # sites are fissionable -- close to OpenMC's default 5 % floor. The initial source only seeds
    # the inactive batches, so a lower acceptance floor is harmless.
    s.source_rejection_fraction = 0.01
    s.temperature = {'method': 'interpolation', 'tolerance': 200.0}
    s.output = {'tallies': False}
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rods', choices=['in', 'out'], default='in')
    ap.add_argument('--particles', type=int, default=20000)
    ap.add_argument('--batches', type=int, default=150)
    ap.add_argument('--inactive', type=int, default=50)
    ap.add_argument('--fuel-T', type=float, default=gfmr.T_FUEL_MEAN_C)
    ap.add_argument('--graphite-T', type=float, default=gfmr.T_GRAPHITE_C)
    ap.add_argument('--enrichment', type=float, default=gfmr.ENRICHMENT)
    ap.add_argument('--out', default='run')
    args = ap.parse_args()

    mats = materials(args.fuel_T, args.graphite_T, args.enrichment)
    geom, _ = geometry(mats, rods_in=(args.rods == 'in'))
    model = openmc.Model(geometry=geom, materials=openmc.Materials(mats.values()),
                         settings=settings(args.particles, args.batches, args.inactive))
    os.makedirs(args.out, exist_ok=True)
    model.export_to_model_xml(os.path.join(args.out, 'model.xml'))
    print('wrote', os.path.join(args.out, 'model.xml'))


if __name__ == '__main__':
    main()
