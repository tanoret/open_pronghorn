# MSTDBTCCorrosionUserObject

`MSTDBTCCorrosionUserObject` evaluates a static Cr/Fe/Ni corrosion endpoint using standard-state
Gibbs functions from an external fluoride/chloride MSTDB-TC pair and the calibrated engineering
closures in `advanced_corrosion_models.json`.

`material_class` denotes the donor/source alloy throughout this endpoint. It is the C++ runtime
mapping of the Python calibration context's `inventory_source_material`; the model does not use it
as a separate cold-leg or deposition-surface material selector.

The raw MSTDB-TC files are not distributed with OpenPronghorn. Both `fluoride_database` and
`chloride_database` are required `FileName` inputs. The wrapper reads the required database edition
and both SHA-256 digests from the advanced parameter database and rejects a version or digest
mismatch. The published parameter set is bound to MSTDB-TC V3.1. `allow_extrapolation` defaults to
false. The advanced artifact also identifies `data/corrosion_database.json` by SHA-256 and stores
the exact base parameters, redox offsets, material densities, and explicit Cr/Fe/Ni transport and
Butler--Volmer properties consumed by the endpoint and its Action bridges. The core compares those
values semantically and rejects a modified or fallback-only base database.

Temperature, exposure, flow, area-to-salt-mass ratio, inventory coupling, material, salt, and redox
are explicit physical inputs. `inventory_coupling_factor` scales dissolved salt inventory only; the
static cold-capture `mass_gain_mg_cm2` closure is independent of it. This distinction belongs to the
frozen calibrated model law; changing it requires a new model revision, refit, and validation. The
wrapper does not select physics from a measurement identifier, source citation, experiment family,
or requested response.

This object evaluates standard-state/Nernst relationships with calibrated activity, transport, and
inventory closures. It does **not** parse ChemSage solution-model declarations, minimize Gibbs
energy, evaluate SUBQ activities, or invoke Thermochimica. Results must therefore be described as an
MSTDB-TC standard-state reduced corrosion model, not as a full thermochemical equilibrium result.

Use [AdvancedCorrosionModelPostprocessor](../postprocessors/AdvancedCorrosionModelPostprocessor.md)
to write any scalar endpoint or species diagnostic.

!syntax parameters /UserObjects/MSTDBTCCorrosionUserObject

!syntax inputs /UserObjects/MSTDBTCCorrosionUserObject

!syntax children /UserObjects/MSTDBTCCorrosionUserObject
