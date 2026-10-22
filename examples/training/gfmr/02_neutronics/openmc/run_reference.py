"""GFMR OpenMC reference calculations (Monte Carlo, continuous energy, ENDF/B-VIII.0).

Cases (each in its own sub-directory of --out):
  heating       nominal with coupled neutron-photon transport: energy deposition salt / graphite
  nominal       rods at 36 cm (SPEC), 650 C fuel / 800 C graphite: k-eff, power map, graphite
                heating, fast (> 0.1 MeV) flux in the graphite, kinetics (IFP beta_eff, Lambda),
                and the multigroup cross sections for the MOOSE diffusion model
  rods_out      all rods withdrawn                          -> worth of the nominal insertion
  rods_in       all rods fully inserted (3.6 m)              -> total rod worth, shutdown margin
  fuel_hot      fuel +100 K (density from the SPEC law)     -> fuel temperature coefficient
  graphite_hot  graphite +100 K                             -> graphite temperature coefficient

Usage: python run_reference.py --out DIR [--particles N] [--batches B] [--inactive I]
                                [--cases nominal rods_out ...] [--threads T]
"""

import argparse
import math
import os
import sys

import numpy as np
import openmc
import openmc.mgxs

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'common'))
import gfmr  # noqa: E402
import build_model as bm  # noqa: E402

# 8-group CASMO structure (eV): thermal groups resolve the graphite thermalization
GROUPS = openmc.mgxs.EnergyGroups('CASMO-8')
N_RADIAL = 35     # mesh tallies: 35 x 35 x 36 (10 cm cells) over the full core
N_AXIAL = 36


def add_tallies(model, cells, mgxs=True):
    tallies = openmc.Tallies()
    H = gfmr.HEIGHT * 100.0
    R = gfmr.R_REFLECTOR * 100.0

    mesh = openmc.RegularMesh(name='core_mesh')
    mesh.lower_left = (-R, -R, 0.0)
    mesh.upper_right = (R, R, H)
    mesh.dimension = (N_RADIAL, N_RADIAL, N_AXIAL)

    # heating (fission fragments + neutron + photon KERMA) on the regular mesh, split by material
    mat_filter = openmc.MaterialFilter([model.materials[0], model.materials[1]])  # fuel, graphite
    t = openmc.Tally(name='heating_mesh')
    t.filters = [openmc.MeshFilter(mesh), mat_filter]
    t.scores = ['heating-local', 'kappa-fission', 'fission']
    tallies.append(t)

    # fast flux (E > 0.1 MeV) and total flux on the mesh, graphite only -> fluence lifetime
    t = openmc.Tally(name='fast_flux_graphite')
    t.filters = [openmc.MeshFilter(mesh), openmc.MaterialFilter([model.materials[1]]),
                 openmc.EnergyFilter([1.0e5, 2.0e7])]
    t.scores = ['flux']
    tallies.append(t)

    # per-channel fission power: one bin per lattice instance of each salt cell (distribcell)
    for c in cells['salt']:
        t = openmc.Tally(name='channel_power_' + c.name)
        t.filters = [openmc.DistribcellFilter(c)]
        t.scores = ['kappa-fission']
        tallies.append(t)

    # global energy deposition split
    t = openmc.Tally(name='global_heating')
    t.filters = [openmc.MaterialFilter([m for m in model.materials
                                        if m.id in model.geometry.get_all_materials()])]
    t.scores = ['heating-local', 'kappa-fission', 'heating']
    tallies.append(t)

    # axial and radial (cylindrical) flux shapes in 8 groups for the diffusion verification
    cmesh = openmc.CylindricalMesh(r_grid=np.linspace(0.0, R, 36), z_grid=np.linspace(0.0, H, 37),
                                   phi_grid=[0.0, 2 * math.pi], name='rz_mesh')
    t = openmc.Tally(name='rz_flux')
    t.filters = [openmc.MeshFilter(cmesh), openmc.EnergyFilter(GROUPS.group_edges)]
    t.scores = ['flux', 'nu-fission', 'kappa-fission']
    tallies.append(t)

    lib = None
    if mgxs:
        lib = openmc.mgxs.Library(model.geometry)
        lib.energy_groups = GROUPS
        lib.num_delayed_groups = 6
        lib.mgxs_types = ['transport', 'nu-transport', 'total', 'absorption', 'nu-fission',
                          'kappa-fission', 'fission', 'nu-scatter matrix', 'chi', 'chi-prompt',
                          'chi-delayed', 'inverse-velocity', 'beta', 'decay-rate',
                          'delayed-nu-fission']
        # One domain per unit-cell region, so that the heterogeneous 3-D diffusion model on the
        # CFD mesh gets its own constants for the salt and the graphite of the inner, outer and
        # control-rod cells, the guide tube, the absorber and the reflector. A cell domain
        # integrates over all lattice instances of that cell (the whole core).
        lib.domain_type = 'cell'
        lib.domains = [c for c in model.geometry.get_all_material_cells().values()]
        lib.correction = 'P0'
        lib.by_nuclide = False
        lib.build_library()
        lib.add_to_tallies(tallies, merge=True)
    model.tallies = tallies
    return lib


def run_case(name, out, args, rods_in=True, insertion=None, dT_fuel=0.0, dT_gr=0.0, mgxs=True,
             ifp=False, photon=False):
    # each case builds a fresh model: reset the auto-generated IDs so that the cases are
    # independent and the tally / filter IDs do not clash
    openmc.reset_auto_ids()
    d = os.path.join(out, name)
    os.makedirs(d, exist_ok=True)
    mats = bm.materials(gfmr.T_FUEL_MEAN_C + dT_fuel, gfmr.T_GRAPHITE_C + dT_gr)
    geom, cells = bm.geometry(mats, rods_in=rods_in, insertion=insertion)
    model = openmc.Model(geometry=geom, materials=openmc.Materials(mats.values()),
                         settings=bm.settings(args.particles, args.batches, args.inactive))
    lib = add_tallies(model, cells, mgxs=mgxs)
    if photon:
        # coupled neutron-photon transport: prompt-gamma heating is deposited where the photons
        # are absorbed (largely in the graphite) instead of at the fission site
        model.settings.photon_transport = True
    if ifp:
        # iterated fission probability: adjoint-weighted beta_eff and generation time
        model.settings.ifp_n_generation = 10
        # total kinetics: one unfiltered tally holding all three IFP scores (StatePoint helper)
        t = openmc.Tally(name='ifp')
        t.scores = ['ifp-time-numerator', 'ifp-beta-numerator', 'ifp-denominator']
        model.tallies.append(t)
        # group-wise beta_eff numerators (divided by the total denominator in postprocess.py)
        t2 = openmc.Tally(name='ifp_groups')
        t2.scores = ['ifp-beta-numerator']
        t2.filters = [openmc.DelayedGroupFilter(list(range(1, 7)))]
        model.tallies.append(t2)
    # MPI ranks scale far better than threads with many tallies (thread contention on the
    # tally arrays): use --ranks R --threads T (R x T cores)
    mpi = ['mpiexec', '-n', str(args.ranks)] if args.ranks > 1 else None
    sp = model.run(cwd=d, threads=args.threads, mpi_args=mpi)
    with openmc.StatePoint(sp) as s:
        k = s.keff
        print('%-14s k-eff = %.5f +/- %.5f' % (name, k.nominal_value, k.std_dev), flush=True)
    if mgxs:
        lib.dump_to_file('mgxs_library', directory=d)
    return k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='reference')
    ap.add_argument('--particles', type=int, default=200000)
    ap.add_argument('--batches', type=int, default=300)
    ap.add_argument('--inactive', type=int, default=100)
    ap.add_argument('--threads', type=int, default=None)
    ap.add_argument('--ranks', type=int, default=1, help='MPI ranks (needs an MPI build of OpenMC)')
    ap.add_argument('--cases', nargs='+',
                    default=['nominal', 'heating', 'rods_out', 'fuel_hot', 'graphite_hot'])
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    spec = dict(nominal=dict(ifp=True), heating=dict(mgxs=False, photon=True), rods_out=dict(rods_in=False),
                rods_in=dict(insertion=gfmr.HEIGHT, mgxs=False),
                fuel_hot=dict(dT_fuel=100.0), graphite_hot=dict(dT_gr=100.0))
    for c in args.cases:
        run_case(c, args.out, args, **spec[c])


if __name__ == '__main__':
    main()
