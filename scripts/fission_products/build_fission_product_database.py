#!/usr/bin/env python3
"""Build data/fission_product_database.json from the ENDF/B-VIII.0 sub-libraries.

Usage:
    python3 build_fission_product_database.py <endf_dir> [output.json]

<endf_dir> must contain the unzipped NNDC sub-libraries
    ENDF-B-VIII.0_decay/    (https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_decay.zip)
    ENDF-B-VIII.0_nfy/      (https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_nfy.zip)
    ENDF-B-VIII.0_neutrons/ (only the n-*.endf files of the nuclides below are read)

What is taken from where:
  * half-lives, decay modes and branching ratios: decay sub-library, MF8/MT457;
  * independent (MT454) and cumulative (MT459) U-235 fission yields at 0.0253 eV ("thermal")
    and 500 keV ("fast"): NFY sub-library;
  * 2200 m/s capture cross sections: the value tabulated by the evaluator in the MF1/MT451
    description of each neutron evaluation (Mughabghab, Atlas of Neutron Resonances), which is
    the value the evaluated MF3/MT102 data reproduce at 0.0253 eV;
  * atomic masses: AWR (MF8 head record) times the neutron mass in u.

Short-lived states that are not tracked (isomers such as Xe-135m, Tc-99m is tracked) are
"lumped": their decay branches are followed through to the next tracked nuclide, so a branching
ratio in the output is the probability that the parent eventually reaches that daughter before
any other tracked nuclide.
"""
import json
import math
import os
import re
import sys

NEUTRON_MASS_U = 1.00866491595
RTYP = {1: "beta-", 2: "ec", 3: "it", 4: "alpha", 5: "n", 6: "sf", 7: "p"}

# Tracked nuclides: (name, Z, symbol, A, isomeric state, chemical group)
TRACKED = [
    ("Kr85", 36, "Kr", 85, 0, "noble_gas"),
    ("Sr90", 38, "Sr", 90, 0, "salt_seeking"),
    ("Y90", 39, "Y", 90, 0, "salt_seeking"),
    ("Zr95", 40, "Zr", 95, 0, "salt_seeking"),
    ("Nb95", 41, "Nb", 95, 0, "noble_metal"),
    ("Mo95", 42, "Mo", 95, 0, "noble_metal"),
    ("Mo99", 42, "Mo", 99, 0, "noble_metal"),
    ("Tc99m", 43, "Tc", 99, 1, "noble_metal"),
    ("Tc99", 43, "Tc", 99, 0, "noble_metal"),
    ("Ru103", 44, "Ru", 103, 0, "noble_metal"),
    ("Rh103", 45, "Rh", 103, 0, "noble_metal"),
    ("Ru106", 44, "Ru", 106, 0, "noble_metal"),
    ("Te132", 52, "Te", 132, 0, "noble_metal"),
    ("I131", 53, "I", 131, 0, "halogen"),
    ("I132", 53, "I", 132, 0, "halogen"),
    ("I133", 53, "I", 133, 0, "halogen"),
    ("Xe133", 54, "Xe", 133, 0, "noble_gas"),
    ("I135", 53, "I", 135, 0, "halogen"),
    ("Xe135", 54, "Xe", 135, 0, "noble_gas"),
    ("Cs135", 55, "Cs", 135, 0, "salt_seeking"),
    ("Cs137", 55, "Cs", 137, 0, "salt_seeking"),
    ("Pm149", 61, "Pm", 149, 0, "salt_seeking"),
    ("Sm149", 62, "Sm", 149, 0, "salt_seeking"),
]

# Capture product (n,gamma) of each tracked nuclide (A+1, ground state)
SYMBOLS = {36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc",
           44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 51: "Sb", 52: "Te", 53: "I", 54: "Xe",
           55: "Cs", 56: "Ba", 57: "La", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu"}


def fnum(s):
    s = s.strip()
    if not s:
        return 0.0
    for i in range(1, len(s)):
        if s[i] in "+-" and s[i - 1] not in "eE":
            return float(s[:i] + "e" + s[i:])
    return float(s)


def records(path, mf, mt):
    out = []
    with open(path) as f:
        for line in f:
            if len(line) < 75:
                continue
            try:
                if int(line[70:72]) == mf and int(line[72:75]) == mt:
                    out.append(line[:66])
            except ValueError:
                pass
    return out


def fields(lines):
    return [l[11 * k:11 * k + 11] for l in lines for k in range(6)]


def name_of(z, a, m):
    return f"{SYMBOLS.get(z, 'Z%d' % z)}{a}" + ("m" if m else "")


def decay_file(root, z, a, m):
    tag = f"{a:03d}" + (f"m{m}" if m else "")
    return os.path.join(root, "ENDF-B-VIII.0_decay", f"dec-{z:03d}_{SYMBOLS[z]}_{tag}.endf")


def read_decay(root, z, a, m):
    path = decay_file(root, z, a, m)
    if not os.path.exists(path):
        return None
    lines = records(path, 8, 457)
    head = lines[0]
    awr = fnum(head[11:22])
    nst = int(head[44:55])
    t12 = fnum(lines[1][0:11])
    if nst == 1:
        return dict(awr=awr, half_life=None, modes=[])
    nc2 = int(lines[1][44:55])
    idx = 2 + (nc2 + 5) // 6
    ndk = int(lines[idx][55:66])
    v = fields(lines[idx + 1: idx + 1 + (6 * ndk + 5) // 6])
    modes = []
    for k in range(ndk):
        rtyp = fnum(v[6 * k])
        rfs = int(fnum(v[6 * k + 1]))
        br = fnum(v[6 * k + 4])
        # only single-step modes (beta-, EC, IT) are needed for the fission products here
        first = int(rtyp)
        dz, da = {1: (1, 0), 2: (-1, 0), 3: (0, 0), 4: (-2, -4), 5: (0, -1), 7: (-1, -1)}.get(
            first, (0, 0))
        modes.append(dict(mode=RTYP.get(first, str(rtyp)), z=z + dz, a=a + da, m=rfs, br=br))
    return dict(awr=awr, half_life=t12, modes=modes)


def read_yields(root):
    path = os.path.join(root, "ENDF-B-VIII.0_nfy", "nfy-092_U_235.endf")
    result = {}
    for mt, label in ((454, "independent"), (459, "cumulative")):
        lines = records(path, 8, mt)
        n_e = int(lines[0][22:33])
        pos = 1
        for _ in range(n_e):
            e = fnum(lines[pos][0:11])
            nfp = int(lines[pos][55:66])
            nl = (4 * nfp + 5) // 6
            v = fields(lines[pos + 1: pos + 1 + nl])
            spectrum = "thermal" if e < 1.0 else ("fast" if e < 1.0e6 else None)
            if spectrum:
                for k in range(nfp):
                    za = int(fnum(v[4 * k]))
                    st = int(fnum(v[4 * k + 1]))
                    result[(za, st, spectrum, label)] = fnum(v[4 * k + 2])
            pos += 1 + nl
    return result


def read_capture(root, z, sym, a, m):
    if m:
        return None
    path = os.path.join(root, "ENDF-B-VIII.0_neutrons", f"n-{z:03d}_{sym}_{a:03d}.endf")
    if not os.path.exists(path):
        return None
    desc = [l[:66] for l in open(path) if l[70:75] == " 1451"]
    vals = []
    for l in desc:
        mm = re.match(r"\s*capture\s+([0-9.+\-Ee]+)", l, re.I)
        if mm:
            vals.append(fnum(mm.group(1)))
    # The 2200 m/s table appears last when an evaluation lists both a resonance-parameter and a
    # pointwise value; they agree for every nuclide used here except Ru-103, whose final (adopted)
    # value is used.
    return vals[-1] if vals else None


def main():
    root = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else "fission_product_database.json"
    tracked = {(z, a, m): n for n, z, s, a, m, g in TRACKED}
    ylds = read_yields(root)

    def lumped_daughters(z, a, m, weight, seen):
        """Follow decays of an untracked state until tracked (or untracked stable) nuclides."""
        d = read_decay(root, z, a, m)
        out = {}
        if d is None or d["half_life"] is None:
            out[name_of(z, a, m)] = weight
            return out
        for mode in d["modes"]:
            key = (mode["z"], mode["a"], mode["m"])
            w = weight * mode["br"]
            if key in tracked:
                out[tracked[key]] = out.get(tracked[key], 0.0) + w
            elif mode["m"] > 0 and key not in seen:
                for k, v in lumped_daughters(*key, w, seen | {key}).items():
                    out[k] = out.get(k, 0.0) + v
            else:
                n = name_of(*key)
                out[n] = out.get(n, 0.0) + w
        return out

    nuclides = {}
    for name, z, sym, a, m, group in TRACKED:
        d = read_decay(root, z, a, m)
        entry = dict(Z=z, A=a, isomeric_state=m, group=group,
                     molar_mass=round(d["awr"] * NEUTRON_MASS_U, 6))
        if d["half_life"] is None:
            entry["half_life"] = None
            entry["decay"] = []
        else:
            entry["half_life"] = d["half_life"]
            daughters = {}
            for mode in d["modes"]:
                key = (mode["z"], mode["a"], mode["m"])
                if key in tracked:
                    daughters[tracked[key]] = daughters.get(tracked[key], 0.0) + mode["br"]
                elif mode["m"] > 0:
                    for k, v in lumped_daughters(*key, mode["br"], {key}).items():
                        daughters[k] = daughters.get(k, 0.0) + v
                else:
                    n = name_of(*key)
                    daughters[n] = daughters.get(n, 0.0) + mode["br"]
            entry["decay"] = [dict(daughter=k, branching=round(v, 8))
                              for k, v in sorted(daughters.items(), key=lambda kv: -kv[1])]
        sigma = read_capture(root, z, sym, a, m)
        entry["capture"] = dict(cross_section=sigma if sigma is not None else 0.0,
                                product=name_of(z, a + 1, 0))
        za = 1000 * z + a
        entry["yields"] = {
            f"U235_{sp}": dict(independent=ylds.get((za, m, sp, "independent"), 0.0),
                               cumulative=ylds.get((za, m, sp, "cumulative"), 0.0))
            for sp in ("thermal", "fast")
        }
        nuclides[name] = entry

    db = dict(
        description="Fission-product chains for the OpenPronghorn FissionProductInventory action.",
        units=dict(half_life="s", cross_section="barn (2200 m/s capture)", molar_mass="g/mol",
                   yields="fraction per fission"),
        provenance=dict(
            decay="ENDF/B-VIII.0 decay sub-library, MF8/MT457 (half-lives, branching ratios).",
            yields="ENDF/B-VIII.0 neutron fission yield sub-library, U-235, MT454 (independent) and "
                   "MT459 (cumulative) at 0.0253 eV (thermal) and 500 keV (fast).",
            capture="2200 m/s capture cross sections tabulated in the MF1/MT451 description of the "
                    "ENDF/B-VIII.0 neutron evaluations (Mughabghab, Atlas of Neutron Resonances).",
            lumping="Untracked isomeric states are lumped: their decay branches are followed to "
                    "the next tracked nuclide (e.g. I135 -> Xe135m -> Xe135).",
            generator="scripts/fission_products/build_fission_product_database.py"),
        groups=dict(
            noble_gas="Insoluble noble gases (Kr, Xe): removed by sparging/off-gas.",
            noble_metal="Noble and semi-noble metals (Nb, Mo, Tc, Ru, Rh, Te): insoluble in the "
                        "salt as metals; plate on metal surfaces and are entrained in the off-gas.",
            halogen="Halogens (I): remain in the salt as iodide.",
            salt_seeking="Salt-seeking fission products (Sr, Y, Zr, Cs, Pm, Sm): stable fluorides "
                         "or chlorides that remain dissolved."),
        nuclides=nuclides,
    )
    with open(out_path, "w") as f:
        json.dump(db, f, indent=2)
        f.write("\n")
    for n, e in nuclides.items():
        hl = e["half_life"]
        print(f"{n:6s} T1/2={'stable' if hl is None else '%.4e s' % hl:>14s} "
              f"sigma={e['capture']['cross_section']:.4g} b "
              f"Ycum={e['yields']['U235_thermal']['cumulative']:.4e} "
              f"decay={[(d['daughter'], d['branching']) for d in e['decay']]}")


if __name__ == "__main__":
    main()
