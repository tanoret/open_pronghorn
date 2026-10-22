# Stage 4: radiolysis and corrosion — NOT DONE

This stage is planned but **not implemented or validated**. Do not use it in the training as a
result.

## Plan

Run the existing OpenPronghorn chemistry actions on the frozen stage-3 flow field, as in case 05
of `examples/training/chemistry` (the MSFR capstone):

- `[MoltenSaltRadiolysis]`: the radiolytic dose taken from the stage-2 power density (about
  0.87 × the fission power in a fuel salt);
- `[CorrosionPlatingFlow]`: element-resolved Hastelloy N corrosion and plating on the guide-tube
  walls (`wall_guide`), with temperature-dependent equilibria from the CFD wall temperature;
- `[FissionProductInventory]`: I-135 / Xe-135 and noble metals, with graphite uptake as a sink.

## Open questions before implementation

- The GFMR specification has no loop, off-gas system or redox-control parameters. These must
  be assumed and stated.
- FUNaK has no U(III)/U(IV) redox data comparable to FLiBe. The redox potential must be
  assumed.
- The only metal in the core is the guide tube (Hastelloy N). Loop corrosion needs a loop model.
