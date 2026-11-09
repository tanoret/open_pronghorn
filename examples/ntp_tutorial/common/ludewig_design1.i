# Ludewig et al. (1996) Design 1 particle-bed element, P/D = 1.75: geometry and materials.
#
# Radii from Table 4.2.4 of Ludewig et al., Prog. Nucl. Energy 30(1) (1996), as used in the
# reference channel model: index 1..5 are the hot-frit inner, hot-frit outer, fuel-bed outer,
# cold-frit outer and plenum outer radii. '_in' is the narrow inlet end (z = 0), '_out' the
# chamber exit (z = L). The reference channel model widens the plenum outer radius at z = L to 33.75 mm
# (Ludewig: 34.75 mm); that value is kept so the two models share one geometry.

L = 0.7294
r1_in = 0.00801
r1_out = 0.02129
r2_in = 0.01101
r2_out = 0.02429
r3_in = 0.01976
r3_out = 0.02975
r4_in = 0.02176
r4_out = 0.03175
r5_in = 0.02376
r5_out = 0.03375

n_elements = 37

# Porous regions (as in the reference channel model)
# cold frit: beryllium, 80 um particles
eps_cf = 0.32
d_cf = 8.0e-5
k_cf = 190.0
# fuel bed: coated fuel particles, 500 um
eps_bed = 0.40
d_bed = 5.0e-4
k_bed = 22.0
kz_bed = 0.65
kr_bed = 0.85
# hot frit: TaC, 120 um
eps_hf = 0.28
d_hf = 1.2e-4
k_hf = 33.0

# Free-flow friction factor of the plenum and chamber annuli
f_ann = 0.02
