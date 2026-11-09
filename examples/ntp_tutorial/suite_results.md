# NTP tutorial: demonstration suite results

Written by `run_suite.py` on 2026-11-09 14:56, executable `../../open_pronghorn-opt`. 129 checks against the READMEs: 129 match, 0 mismatch or missing, 0 not run; 0 job(s) failed, 0 skipped.

Tolerance: a documented number matches when the fresh value rounds to it; round-off quantities are upper bounds; run times within a factor of three. The CSV `suite_results.csv` holds every extracted figure (SI units) including the unchecked ones.

## Jobs

| job | group | status | ranks | wall [s] | command (from the stage directory) |
|---|---|---|---|---|---|
| 01_hydrogen | 01 | ok | 1 | 1.50548 | `cd 01_closures && python3 check_hydrogen.py --exe open_pronghorn-opt` |
| 01_validate | 01 | ok | 1 | 66.058 | `cd 01_closures && python3 validate.py --exe open_pronghorn-opt` |
| 02_mesh | 02 | ok | 1 | 3.00466 | `cd 02_element && open_pronghorn-opt -i mesh.i --mesh-only --allow-unused --color off` |
| 02_plot_mesh | 02 | ok | 1 | 1.50045 | `cd 02_element && python3 plot_mesh.py` |
| 02_newton_cold | 02 | ok | 1 | 6.00905 | `cd 02_element && open_pronghorn-opt -i ../03_ludewig_anchor/anchor.i ../03_ludewig_anchor/continuation.i Executioner/end_time=0 Outputs/file_base=newton_cold --color off` |
| 02_newton_5pc | 02 | ok | 1 | 3.5048 | `cd 02_element && open_pronghorn-opt -i ../03_ludewig_anchor/anchor.i ../03_ludewig_anchor/continuation.i Executioner/end_time=0.1 Outputs/file_base=newton_5pc --color off` |
| 02_simple_cold | 02 | ok | 1 | 7.00567 | `cd 02_element && open_pronghorn-opt -i element_simple_rz.i Outputs/file_base=simple_cold --color off` |
| 02_simple_5pc | 02 | ok | 1 | 48.0264 | `cd 02_element && open_pronghorn-opt -i element_simple_rz.i s_power=0.05 simple_fixed_point=40 --color off` |
| 02_compare_cold_flow | 02 | ok | 1 | 1.50043 | `cd 02_element && python3 compare_cold_flow.py` |
| 03_anchor | 03 | ok | 4 | 7.50183 | `cd 03_ludewig_anchor && mpirun --bind-to none -n 4 open_pronghorn-opt -i anchor.i continuation.i orifice_implicit.i Outputs/file_base=anchor_out --color off` |
| 03_fold_uniform | 03 | ok (turning point) | 1 | 46.0208 | `cd 03_ludewig_anchor && open_pronghorn-opt -i anchor.i continuation.i s_design=0.3 t_design=0.6 Executioner/dtmin=1e-6 Executioner/error_on_dtmin=false --color off` — stopped at the turning point by Executioner/dtmin=1e-6 Executioner/error_on_dtmin=false |
| 03_fold_chi30 | 03 | ok (turning point) | 1 | 62.5285 | `cd 03_ludewig_anchor && open_pronghorn-opt -i anchor.i continuation.i s_design=0.3 t_design=0.6 cold_frit_multiplier=30 Outputs/file_base=fold_chi30_out Executioner/dtmin=1e-6 Executioner/error_on_dtmin=false --color off` — uniform throttling chi = 30 (README: fold at 3.3 %) |
| 03_fold_chi100 | 03 | ok (turning point) | 1 | 58.528 | `cd 03_ludewig_anchor && open_pronghorn-opt -i anchor.i continuation.i s_design=0.3 t_design=0.6 cold_frit_multiplier=100 Outputs/file_base=fold_chi100_out Executioner/dtmin=1e-6 Executioner/error_on_dtmin=false --color off` — uniform throttling chi = 100 (README: fold at 2.0 %) |
| 03_anchor_nz24 | 03 | ok | 4 | 11.0029 | `cd 03_ludewig_anchor && mpirun --bind-to none -n 4 open_pronghorn-opt -i anchor.i continuation.i orifice_implicit.i nz=24 log_chi_order=TWENTYFOURTH Outputs/file_base=anchor_nz24_out --color off` — mesh variant of the README; no documented values |
| 04_simple_flow | 04 | ok | 1 | 118.039 | `cd 04_navier_stokes && open_pronghorn-opt -i anchor_simple.i orifice_flow.i orifice_verbose=false Outputs/file_base=anchor_simple_out --color off` — (c) SIMPLE, orifice designed on the slice flows |
| 04_pimple_flow | 04 | ok | 1 | 181.554 | `cd 04_navier_stokes && open_pronghorn-opt -i anchor_pimple.i orifice_flow.i orifice_verbose=false Outputs/file_base=anchor_pimple_out --color off` — (c) PIMPLE, orifice designed on the slice flows |
| 04_extract_orifice | 04 | ok | 1 | 0.500187 | `cd 04_navier_stokes && python3 extract_orifice.py ../03_ludewig_anchor/anchor_out.csv` |
| 04_pimple_newton_orifice | 04 | ok (runaway, as documented) | 1 | 265.047 | `cd 04_navier_stokes && open_pronghorn-opt -i anchor_pimple.i Outputs/file_base=anchor_pimple_newton_out --color off` — (a) finding: the Newton orifice starves the middle slices (README) |
| 04_pimple_lagged | 04 | ok (runaway, as documented) | 1 | 404.088 | `cd 04_navier_stokes && open_pronghorn-opt -i anchor_pimple.i orifice_lagged.i s_final=0.5 t_ramp=30 t_end=60 pimple_steady_state_detection=false Outputs/file_base=orifice_lagged_out --color off` — (b) finding: lagged temperature design along the ramp to 50 % |
| 05_mesh | 05 | ok | 1 | 2.00046 | `cd 05_element_3d && open_pronghorn-opt -i mesh_3d.i --mesh-only --allow-unused --color off` |
| 05_anchor_3d | 05 | ok | 12 | 37.5072 | `cd 05_element_3d && mpirun --bind-to none -n 12 open_pronghorn-opt -i anchor_3d.i ../03_ludewig_anchor/continuation.i ../03_ludewig_anchor/orifice_implicit.i tolerance_3d.i Outputs/file_base=anchor_3d_out --color off` |
| 05_vs_03 | 05 | ok | 1 | 0 | `(comparison of outputs)` — 3-D vs axisymmetric anchor |
| 05_tilt_005 | long | ok (turning point) | 8 | 133.06 | `cd 05_element_3d && mpirun --bind-to none -n 8 open_pronghorn-opt -i anchor_3d.i ../03_ludewig_anchor/continuation.i ../03_ludewig_anchor/orifice_implicit.i tolerance_3d.i tilt=0.05 Outputs/file_base=anchor_3d_tilt_out Executioner/dtmin=1e-6 Executioner/error_on_dtmin=false --color off` — stopped at the turning point by Executioner/dtmin=1e-6 Executioner/error_on_dtmin=false |
| 05_tilt_015 | long | ok (turning point) | 8 | 107.528 | `cd 05_element_3d && mpirun --bind-to none -n 8 open_pronghorn-opt -i anchor_3d.i ../03_ludewig_anchor/continuation.i ../03_ludewig_anchor/orifice_implicit.i tolerance_3d.i tilt=0.15 Outputs/file_base=anchor_3d_tilt015_out Executioner/dtmin=1e-6 Executioner/error_on_dtmin=false --color off` — stopped at the turning point by Executioner/dtmin=1e-6 Executioner/error_on_dtmin=false |
| 05_tilt_001 | long | ok (turning point) | 8 | 164.751 | `cd 05_element_3d && mpirun --bind-to none -n 8 open_pronghorn-opt -i anchor_3d.i ../03_ludewig_anchor/continuation.i ../03_ludewig_anchor/orifice_implicit.i tolerance_3d.i tilt=0.01 Outputs/file_base=anchor_3d_tilt001_out Executioner/dtmin=1e-6 Executioner/error_on_dtmin=false --color off` — stopped at the turning point by Executioner/dtmin=1e-6 Executioner/error_on_dtmin=false |
| 05_simple_cold_3d | 05 | ok | 8 | 15.503 | `cd 05_element_3d && mpirun --bind-to none -n 8 open_pronghorn-opt -i element_simple_3d.i Outputs/file_base=element_simple_3d_out --color off` — cold flow from the plug-flow start |
| 05_simple_3d | long | ok | 8 | 301.771 | `cd 05_element_3d && mpirun --bind-to none -n 8 open_pronghorn-opt -i anchor_simple_3d.i ../04_navier_stokes/orifice_flow.i orifice_sectors=8 orifice_verbose=false Outputs/file_base=anchor_simple_3d_out --color off` — SIMPLE, column orifice (12 slices x 8 sectors) |
| 05_pimple_3d | long | ok | 8 | 488.603 | `cd 05_element_3d && mpirun --bind-to none -n 8 open_pronghorn-opt -i anchor_pimple_3d.i ../04_navier_stokes/orifice_flow.i orifice_sectors=8 orifice_verbose=false Outputs/file_base=anchor_pimple_3d_out --color off` — PIMPLE from rest, column orifice |
| 05_pimple_3d_rings | long | ok | 8 | 184.047 | `cd 05_element_3d && mpirun --bind-to none -n 8 open_pronghorn-opt -i anchor_pimple_3d.i ../04_navier_stokes/orifice_flow.i t_end=10 pimple_steady_state_detection=false orifice_verbose=false Outputs/file_base=anchor_pimple_3d_rings_out --color off` — ring orifice: the azimuthal flow excursion |
| 05_plot_3d | long | ok | 1 | 1.50037 | `cd 05_element_3d && python3 plot_3d.py 'anchor_3d_out.e:no tilt' 'anchor_3d_tilt_out.e:tilt 0.05'` |

## Figures of merit vs the READMEs

| job | quantity | this run | documented | tolerance | status | source |
|---|---|---|---|---|---|---|
| 01_hydrogen | dh_30_3000 | 4.750e+07 | 4.749625e7 J/kg | ±5 J/kg | ok | 01_closures/README.md |
| 01_hydrogen | mdot_bundle | 21.0543 | 21.054295 kg/s | ±5e-07 kg/s | ok | 01_closures/README.md |
| 01_hydrogen | mdot_element | 0.569035 | 0.569035 kg/s | ±5e-07 kg/s | ok | 01_closures/README.md |
| 01_hydrogen | props_max_rel_diff | 6.990e-14 | 7e-14 | |x| <= 1e-12 | ok | 01_closures/README.md |
| 01_hydrogen | props_max_rel_diff | 6.990e-14 | 7e-14 | |x| <= 1e-12 | ok | README.md |
| 01_validate | data_max_err_pc | 13.73 | 13.73 % | ±0.005 % | ok | 01_closures/README.md |
| 01_validate | data_max_err_pc | 13.73 | 13.7 % | ±0.05 % | ok | README.md |
| 01_validate | data_mean_err_pc | 5.32 | 5.32 % | ±0.005 % | ok | 01_closures/README.md |
| 01_validate | data_mean_err_pc | 5.32 | 5.3 % | ±0.05 % | ok | README.md |
| 01_validate | newton_vs_closed_form | 8.240e-12 | 8e-12 | |x| <= 1e-10 | ok | 01_closures/README.md |
| 01_validate | points | 46 | 46 | ±0.5 | ok | 01_closures/README.md |
| 01_validate | simple_vs_closed_form | 1.770e-13 | 2e-13 | |x| <= 1e-11 | ok | 01_closures/README.md |
| 02_newton_cold | dp_in | 16.2084 | 16.2 kPa | ±0.05 kPa | ok | 02_element/README.md |
| 02_newton_cold | mass_error | -9.401e-07 | -9e-7 | ±5e-08 | ok | 02_element/README.md |
| 02_newton_5pc | T_fuel_max | 324.712 | 325 K | ±0.5 K | ok | 02_element/README.md |
| 02_newton_5pc | T_mixed | 218.922 | 218.9 K | ±0.05 K | ok | 02_element/README.md |
| 02_newton_5pc | energy_error | 8.442e-08 | 8e-8 | ±5e-09 | ok | 02_element/README.md |
| 02_simple_cold | dp_in | 11.5538 | 11.6 kPa | ±0.05 kPa | ok | 02_element/README.md |
| 02_simple_cold | mass_error | 4.706e-13 | 5e-13 | |x| <= 1e-11 | ok | 02_element/README.md |
| 02_simple_5pc | T_fuel_max | 1085.42 | 1085 K | ±0.5 K | ok | 02_element/README.md |
| 02_simple_5pc | T_mixed | 218.948 | 218.9 K | ±0.05 K | ok | 02_element/README.md |
| 02_simple_5pc | energy_error | -5.217e-08 | -5e-8 | ±5e-09 | ok | 02_element/README.md |
| 03_anchor | T_fuel_max | 3193.17 | 3193 K | ±0.5 K | ok | 03_ludewig_anchor/README.md |
| 03_anchor | T_fuel_max | 3193.17 | 3193 K | ±0.5 K | ok | README.md |
| 03_anchor | T_fuel_max | 3193.17 | 3193.2 K | ±0.05 K | ok | 05_element_3d/README.md |
| 03_anchor | T_mixed | 2989.42 | 2989.4 K | ±0.05 K | ok | 03_ludewig_anchor/README.md |
| 03_anchor | T_mixed | 2989.42 | 2989.4 K | ±0.05 K | ok | README.md |
| 03_anchor | T_mixed | 2989.42 | 2989.4 K | ±0.05 K | ok | 05_element_3d/README.md |
| 03_anchor | chi_exit | 20.6305 | 20.6 | ±0.05 | ok | 03_ludewig_anchor/README.md |
| 03_anchor | chi_exit | 20.6305 | 20.63 | ±0.005 | ok | 05_element_3d/README.md |
| 03_anchor | chi_inlet | 0.147981 | 0.15 | ±0.005 | ok | 03_ludewig_anchor/README.md |
| 03_anchor | energy_error | 3.189e-07 | 3e-7 | ±5e-08 | ok | 03_ludewig_anchor/README.md |
| 03_anchor | energy_error | 3.189e-07 | 3e-7 | ±5e-08 | ok | README.md |
| 03_anchor | energy_error | 3.189e-07 | 3e-7 | ±5e-08 | ok | 05_element_3d/README.md |
| 03_anchor | mass_error | 7.082e-14 | 1e-13 | |x| <= 1e-12 | ok | 03_ludewig_anchor/README.md |
| 03_anchor | mass_error | 7.082e-14 | 1e-13 | |x| <= 1e-12 | ok | README.md |
| 03_anchor | p_in | 7.84118 | 7.84 MPa | ±0.005 MPa | ok | 03_ludewig_anchor/README.md |
| 03_anchor | p_in | 7.84118 | 7.84 MPa | ±0.005 MPa | ok | README.md |
| 03_anchor | p_in | 7.84118 | 7.841 MPa | ±0.0005 MPa | ok | 05_element_3d/README.md |
| 03_anchor | power | 26.9164 | 26.92 MW | ±0.005 MW | ok | 03_ludewig_anchor/README.md |
| 03_anchor | run_time | 7.50183 | 10 s | ×/÷3 | ok | 03_ludewig_anchor/README.md ('about 10 s, 4 ranks') |
| 03_anchor | s_power | 1 | 1 | ±0.5 | ok | 03_ludewig_anchor/README.md |
| 03_fold_uniform | s_max | 0.0853931 | 0.085 | ±0.0005 | ok | 03_ludewig_anchor/README.md (8.5 %) |
| 03_fold_uniform | s_max | 0.0853931 | 0.085 | ±0.0005 | ok | README.md (8.5 %) |
| 03_fold_chi30 | s_max | 0.032981 | 0.033 | ±0.0005 | ok | 03_ludewig_anchor/README.md (3.3 % at chi = 30) |
| 03_fold_chi100 | s_max | 0.0196715 | 0.020 | ±0.0005 | ok | 03_ludewig_anchor/README.md (2.0 % at chi = 100) |
| 03_anchor_nz24 | s_power | 1 | 1 | ±0.5 | ok | 03_ludewig_anchor/README.md (reaches the design power) |
| 04_simple_flow | T_fuel_max | 3167.38 | 3167 K | ±0.5 K | ok | 04_navier_stokes/README.md |
| 04_simple_flow | T_fuel_max | 3167.38 | 3167 K | ±0.5 K | ok | README.md |
| 04_simple_flow | T_mixed | 2989.81 | 2989.81 K | ±0.005 K | ok | 04_navier_stokes/README.md |
| 04_simple_flow | T_mixed | 2989.81 | 2989.8 K | ±0.05 K | ok | README.md |
| 04_simple_flow | T_slice_spread | 0.0011195 | 1.1e-3 | ±5e-05 | ok | 04_navier_stokes/README.md |
| 04_simple_flow | chi_max | 46.2682 | 46.3 | ±0.05 | ok | 04_navier_stokes/README.md |
| 04_simple_flow | chi_min | 0.406679 | 0.41 | ±0.005 | ok | 04_navier_stokes/README.md |
| 04_simple_flow | energy_error | -5.393e-09 | -5e-9 | |x| <= 1e-08 | ok | 04_navier_stokes/README.md |
| 04_simple_flow | mass_error | -2.771e-14 | 3e-14 | |x| <= 1e-12 | ok | 04_navier_stokes/README.md |
| 04_simple_flow | p_in | 8.87294 | 8.873 MPa | ±0.0005 MPa | ok | 04_navier_stokes/README.md |
| 04_simple_flow | run_time | 118.039 | 112 s | ×/÷3 | ok | 04_navier_stokes/README.md (1 rank) |
| 04_simple_flow | s_power | 1 | 1 | ±0.5 | ok | 04_navier_stokes/README.md |
| 04_pimple_flow | T_fuel_max | 3167.21 | 3167 K | ±0.5 K | ok | 04_navier_stokes/README.md |
| 04_pimple_flow | T_mixed | 2989.81 | 2989.81 K | ±0.005 K | ok | 04_navier_stokes/README.md |
| 04_pimple_flow | T_slice_spread | 0.00106295 | 1.1e-3 | ±5e-05 | ok | 04_navier_stokes/README.md |
| 04_pimple_flow | chi_max | 46.1399 | 46.1 | ±0.05 | ok | 04_navier_stokes/README.md |
| 04_pimple_flow | chi_min | 0.413446 | 0.41 | ±0.005 | ok | 04_navier_stokes/README.md |
| 04_pimple_flow | energy_error | -7.384e-11 | -7e-11 | |x| <= 1e-08 | ok | 04_navier_stokes/README.md |
| 04_pimple_flow | mass_error | -8.273e-14 | 8e-14 | |x| <= 1e-12 | ok | 04_navier_stokes/README.md |
| 04_pimple_flow | p_in | 8.87219 | 8.872 MPa | ±0.0005 MPa | ok | 04_navier_stokes/README.md |
| 04_pimple_flow | run_time | 181.554 | 165 s | ×/÷3 | ok | 04_navier_stokes/README.md (1 rank) |
| 04_pimple_flow | s_power | 1 | 1 | ±0.5 | ok | 04_navier_stokes/README.md |
| 04_extract_orifice | chi_max | 20.63 | 20.6 | ±0.05 | ok | 03_ludewig_anchor/README.md |
| 04_extract_orifice | chi_min | 0.148 | 0.15 | ±0.005 | ok | 03_ludewig_anchor/README.md |
| 04_pimple_newton_orifice | T_fuel_peak | 3.202e+05 | > 1e4 K | >= 10000 K | ok | 04_navier_stokes/README.md |
| 04_pimple_newton_orifice | T_mixed | 2989.78 | 2990 K | ±0.5 K | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | T_fuel_max@10% | 413.801 | 414 K | ±0.5 K | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | T_fuel_max@25% | 935.343 | 935 K | ±0.5 K | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | T_fuel_max@50% | 1739.23 | 1739 K | ±0.5 K | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | T_fuel_peak | 3.396e+05 | > 1e4 K | >= 10000 K | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | T_mixed@10% | 383.662 | 384 K | ±0.5 K | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | T_mixed@25% | 872.071 | 872 K | ±0.5 K | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | T_mixed@50% | 1634.19 | 1634 K | ±0.5 K | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | chi_max@10% | 93.8667 | 94 | ±0.5 | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | chi_max@50% | 77.0114 | 77 | ±0.5 | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | chi_min@10% | 17.8626 | 17.9 | ±0.05 | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | chi_min@50% | 24.2978 | 24.3 | ±0.05 | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | p_in@10% | 7.25929 | 7.26 MPa | ±0.005 MPa | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | p_in@25% | 7.6565 | 7.66 MPa | ±0.005 MPa | ok | 04_navier_stokes/README.md |
| 04_pimple_lagged | p_in@50% | 8.34438 | 8.34 MPa | ±0.005 MPa | ok | 04_navier_stokes/README.md |
| 05_anchor_3d | T_fuel_max | 3197.09 | 3197.1 K | ±0.05 K | ok | 05_element_3d/README.md |
| 05_anchor_3d | T_mixed | 2989.42 | 2989.4 K | ±0.05 K | ok | 05_element_3d/README.md |
| 05_anchor_3d | chi_exit | 20.6266 | 20.63 | ±0.005 | ok | 05_element_3d/README.md |
| 05_anchor_3d | energy_error | 3.170e-07 | 3e-7 | ±5e-08 | ok | 05_element_3d/README.md |
| 05_anchor_3d | faceting_pc | 2.53582 | 2.5 % | ±0.05 % | ok | 05_element_3d/README.md (faceted bed 2.5 % smaller, n_theta = 8) |
| 05_anchor_3d | fuel_rise_ratio | 1.07026 | 1.07 | ±0.005 | ok | 05_element_3d/README.md (no tilt) |
| 05_anchor_3d | p_in | 7.83611 | 7.836 MPa | ±0.0005 MPa | ok | 05_element_3d/README.md |
| 05_anchor_3d | run_time | 37.5072 | 40 s | ×/÷3 | ok | 05_element_3d/README.md ('~40 s', 12 ranks) |
| 05_anchor_3d | s_power | 1 | 1 | ±0.5 | ok | 05_element_3d/README.md (tilt 0: 1.00) |
| 05_vs_03 | fuel_diff_pc | 0.122898 | 0.12 % | ±0.005 % | ok | 05_element_3d/README.md |
| 05_vs_03 | fuel_diff_pc | 0.122898 | 0.12 % | ±0.005 % | ok | README.md |
| 05_vs_03 | p_in_diff_pc | 0.0646636 | 0.06 % | ±0.005 % | ok | 05_element_3d/README.md |
| 05_tilt_005 | fuel_rise_ratio | 2.16306 | 2.2 | ±0.05 | ok | 05_element_3d/README.md |
| 05_tilt_005 | gas_rise_max | 0.468845 | 0.47 | ±0.005 | ok | 05_element_3d/README.md (+47 %) |
| 05_tilt_005 | gas_rise_max | 0.468845 | 0.47 | ±0.005 | ok | README.md (+47 %) |
| 05_tilt_005 | gas_rise_min | -0.24107 | -0.24 | ±0.005 | ok | 05_element_3d/README.md (-24 %) |
| 05_tilt_005 | gas_rise_min | -0.24107 | -0.24 | ±0.005 | ok | README.md (-24 %) |
| 05_tilt_005 | s_max | 0.178166 | 0.178 | ±0.0005 | ok | 05_element_3d/README.md |
| 05_tilt_005 | s_max | 0.178166 | 0.18 | ±0.005 | ok | README.md (18 %) |
| 05_tilt_015 | s_max | 0.0497028 | 0.050 | ±0.0005 | ok | 05_element_3d/README.md |
| 05_tilt_001 | s_max | 0.513497 | 0.513 | ±0.0005 | ok | 05_element_3d/README.md |
| 05_simple_cold_3d | dp_in | 11.9672 | 11.97 kPa | ±0.005 kPa | ok | 05_element_3d/README.md |
| 05_simple_cold_3d | mass_error | 1.619e-10 | 1e-10 | |x| <= 1e-09 | ok | 05_element_3d/README.md |
| 05_simple_3d | T_fuel_max | 3170.25 | 3170 K | ±0.5 K | ok | 05_element_3d/README.md |
| 05_simple_3d | T_mixed | 2989.42 | 2989.42 K | ±0.005 K | ok | 05_element_3d/README.md |
| 05_simple_3d | T_slice_spread | 0.00109966 | 1.1e-3 | ±5e-05 | ok | 05_element_3d/README.md |
| 05_simple_3d | asymmetry | 1.672e-07 | 2e-7 | |x| <= 1e-06 | ok | 05_element_3d/README.md |
| 05_simple_3d | energy_error | 5.810e-09 | 1e-8 | |x| <= 1e-07 | ok | 05_element_3d/README.md |
| 05_simple_3d | mass_error | 6.096e-12 | 6e-12 | |x| <= 1e-10 | ok | 05_element_3d/README.md |
| 05_simple_3d | p_in | 8.90018 | 8.900 MPa | ±0.0005 MPa | ok | 05_element_3d/README.md |
| 05_simple_3d | run_time | 301.771 | 300 s | ×/÷3 | ok | 05_element_3d/README.md (8 ranks) |
| 05_simple_3d | s_power | 1 | 1 | ±0.5 | ok | 05_element_3d/README.md |
| 05_pimple_3d | T_fuel_max | 3170.09 | 3170 K | ±0.5 K | ok | 05_element_3d/README.md |
| 05_pimple_3d | T_mixed | 2989.42 | 2989.42 K | ±0.005 K | ok | 05_element_3d/README.md |
| 05_pimple_3d | T_slice_spread | 0.00104661 | 1.0e-3 | ±5e-05 | ok | 05_element_3d/README.md |
| 05_pimple_3d | asymmetry | 1.668e-07 | 2e-7 | |x| <= 1e-06 | ok | 05_element_3d/README.md |
| 05_pimple_3d | energy_error | -9.871e-11 | 1e-8 | |x| <= 1e-07 | ok | 05_element_3d/README.md |
| 05_pimple_3d | p_in | 8.89936 | 8.899 MPa | ±0.0005 MPa | ok | 05_element_3d/README.md |
| 05_pimple_3d | run_time | 488.603 | 470 s | ×/÷3 | ok | 05_element_3d/README.md (8 ranks) |
| 05_pimple_3d | s_power | 1 | 1 | ±0.5 | ok | 05_element_3d/README.md |
| 05_pimple_3d_rings | T_fuel_peak | 21984 | > 1e4 K | >= 10000 K | ok | 05_element_3d/README.md |
| 05_pimple_3d_rings | asymmetry | 7.3646 | > 1 | >= 1 | ok | 05_element_3d/README.md |
