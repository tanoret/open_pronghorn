# Case 08: where do the noble metals go? Plating, off-gas and decay on the walls

**Module 2 (fission products) | ~20 min | 1 cell, serial, 18 s**

## Objective

Split a chemical group between two removal pathways, and see why the removed inventories must keep
decaying: their daughters can belong to a different chemical group.

## Physics

The noble and semi-noble fission products (Mo, Tc, Ru, Rh, Te, Nb) do not form stable fluorides in
a reducing salt. They exist as metal atoms or colloids that **plate** on metal surfaces (heat
exchanger, piping) and graphite, or are **entrained into the off-gas** by the sparge bubbles. At
MSRE only a few percent of the Mo-99, Ru-103 and Te-132 produced were found in salt samples
(Compere et al., ORNL-4865, 1975).

Two first-order pathways, plated ($r_p$) and off-gas ($r_o$):

$$\frac{dN}{dt} = \frac{YF}{N_A} + \text{feeds} - (\lambda + r_p + r_o)N, \qquad
  \frac{dN_{plated}}{dt} = r_p N - \lambda N_{plated} + \ldots$$

- Mo-99 (66 h) decays to Tc-99m (6 h) and then Tc-99. Plated Mo-99 produces *plated* Tc-99m.
- Te-132 (3.2 d) decays to I-132 (2.3 h). Plated tellurium produces **iodine on the heat-exchanger
  walls**, even though iodine itself is salt-seeking.

## The input, step by step

[noble_metal_partition.i](noble_metal_partition.i):

- `removal_groups = 'noble_metal noble_metal'`, `removal_rates = '${r_plate} ${r_offgas}'` and
  `removal_destinations = 'plated offgas'`: one group with two pathways.
- `r_plate = 5e-4` 1/s corresponds to $k_m A/V$ with $k_m = 10^{-5}$ m/s on 50 m^2 of heat
  exchanger per m^3 of salt. Case 15 computes $k_m$ from the turbulent flow.
- Run with `verbose = true` and look at `I132_plated` and `I132_offgas`. The action created them
  because I-132 is a daughter of a removed nuclide.

## Run

```bash
open_pronghorn-opt -i noble_metal_partition.i
python plot.py
```

## Expected results

![](../figures/case08_noble_metal_partition.png)

- Only 0.5% of the Mo-99 inventory is in the salt. 83% is plated and 17% is in the off-gas,
  the 5:1 rate ratio.
- `Mo99_salt` matches $YF/(N_A(\lambda + r_p + r_o))$.
- 83% of the I-132 is on the walls.

## Exercises

1. Reduce plating 10x (`r_plate=5e-5`). How do the salt fraction and the plated fraction change?
2. Why does the salt fraction of Ru-103 (0.08%) come out smaller than that of Mo-99 (0.5%)?
3. Why does it matter for the source term that I-132 is born on the walls?

<details><summary>Solutions</summary>

1. The salt fraction of Mo-99 rises from 0.5% to 1.9%, and its plated share falls from 83% to 33%
   (I-132 likewise). Fast removal dominates either way: the salt holds little.
2. Ru-103 lives longer (39 d versus 66 h). Its removed inventory keeps accumulating relative to
   the salt inventory, which saturates at $YF/(N_A(\lambda + r))$.
3. Iodine is volatile and radiologically important. Plated tellurium is a source of iodine
   *outside* the salt, in the heat exchanger, where it is not retained by the salt chemistry.

</details>
