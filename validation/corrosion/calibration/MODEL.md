# Scientific and numerical model

## Scope

This is an effective engineering correlation for salt-side alloy dissolution, mass transfer,
deposition/plating, attack depth, and selected inventory responses. It is not a fully resolved
electrochemical model. Most source experiments do not report electrode potentials, dissolved-metal
activities, wetted areas, velocities, or loop salt volumes consistently. Those unresolved effects
are represented by fitted class offsets and closure coefficients.

The model should therefore be used for source-term screening and sensitivity analysis within the
material, salt, temperature, flow, and redox classes represented by the workbook. Extrapolation
outside that domain should carry an explicit model-form uncertainty.

Chromium solid diffusivity is treated separately from the fitted corrosion/plating coefficients.
Material-specific Arrhenius properties are read from `data/corrosion_database.json` under
`solid_diffusivities`. Hastelloy N uses the DeVan chromium tracer-diffusion correlation.
`generic_metal` retains the historical diffusion correlation as an explicitly labeled legacy
fallback.

## Corrosion branch

The kinetic corrosion rate is

```text
ln(r_kin) = ln(r0)
          + Ea_corr/R * (1/T_ref - 1/T)
          + material_offset
          + salt_offset
          + redox_offset
          + position_offset
          + gamma_flow_corr * ln(flow_factor)
          + theta_dT_corr * ln(1 + DeltaT/100).
```

The transport-capacity branch has its own intercept, activation energy, flow exponent, and
thermal-gradient coefficient. The reported uniform rate is their harmonic combination:

```text
r_corr = 1 / (1/r_kin + 1/r_mt).
```

This form approaches the smaller branch smoothly and prevents high electrochemical drive from
creating an unbounded dissolution rate.

## Effective Butler-Volmer interpretation

The redox offset is additive in log-rate space. It can be written as an equivalent Butler-Volmer
overpotential,

```text
redox_offset = alpha*n*F*eta_effective/(R*T).
```

`eta_effective` must not be interpreted as a measured electrode potential. It absorbs unresolved
activity, speciation, and electrochemical-potential effects shared by a named redox class.

## Deposition and auxiliary responses

The deposition branch contains its own intercept, temperature and flow terms, salt offsets,
cold-leg enhancement, and surface-affinity terms. Its redox contribution has the opposite sense
and a reduced empirical magnitude. Additional fitted closures map the common
corrosion/deposition drive to:

- uniform and intergranular attack depth;
- areal mass loss and gain;
- salt chromium and iron inventory changes;
- noble-metal off-gas fraction; and
- tellurium solubility and redox threshold.

Chromium solid-state diffusivity is not one of these fitted closures in the current model. It is an
independent material property supplied from `solid_diffusivities`.

Agreement in one fitted branch does not establish a mechanistic validation of every underlying
process.

## Chromium solid diffusivity

For a material with chromium Arrhenius properties,

```text
D_Cr(T) = D0 * exp[-Q/(R*T)].
```

The production C++ model and the Python calibration model read the same material-property entries
from `data/corrosion_database.json`.

For exact `hastelloy_n`, the chromium diffusivity uses the DeVan tracer-derived correlation. The
temperature range stored with that entry is provenance metadata for the measured range; the model
does not clamp, reject, or switch correlations outside that range.

If no material-specific chromium entry is present, the model uses the explicitly labeled
`generic_metal` legacy fallback. That fallback reproduces the historical 61-parameter diffusion
correlation but is not presented as a source-specific alloy property.

M-005, the ORNL-TM-3063 chromium-diffusion value at 650 C, is retained as a quantitative
`validation_only` target. It is not used to fit the 59 corrosion/plating coefficients. In the
current framework it checks the independently supplied Hastelloy N diffusivity against the
historical surveillance-derived value.

## Relation-aware residuals

The calibration does not convert every literature statement into an artificial exact value:

- Direct and range targets use weighted log residuals.
- Upper and lower bounds use one-sided residuals that are zero when the inequality is satisfied.
- The deposition-ranking row contributes two ordering residuals.
- `validation_only` targets are predicted and scored but do not enter the fit objective.
- Input-only and excluded auxiliary rows remain in the audit tables but do not enter the objective.

For a direct target,

```text
r_i = sqrt(quality_weight_i)
      * ln(prediction_i / target_i)
      / sigma_ln_i.
```

The optimizer minimizes a robust `soft_l1` loss over the data residuals and 59 Gaussian-style prior
residuals:

```text
r_prior,j = prior_strength * (parameter_j - prior_j) / prior_sigma_j.
```

The 37 active calibration rows produce 38 data residuals because the single deposition-ranking row
encodes two inequalities. The default `prior_strength` is `0.35`. Bounds, initial values, prior
centers, and prior widths are defined next to every fitted parameter in
`src/msr_corrosion_bv/model.py`.

## Optimization

The fit uses `scipy.optimize.least_squares` with:

- bounded parameters;
- `loss="soft_l1"`;
- `f_scale=1`;
- Jacobian-based variable scaling;
- a default maximum of 4000 function evaluations; and
- deterministic initial values.

The frozen run terminated successfully after 21 function evaluations. Because the fit is
deterministic, no random seed is involved.

## Interpretation of reported performance

The current fit uses 43 measurement rows, 37 active calibration rows, and 59 fitted
corrosion/plating parameters plus priors. M-005 is withheld from fitting as one validation-only
quantitative target.

Current quantitative summary metrics are:

| Metric | Value |
|---|---:|
| Direct/range calibration targets | 28 |
| Median direct/range factor error | 1.1166 |
| Direct/range targets within factor 2 | 92.9% |
| Direct/range targets within factor 5 | 100% |
| Validation-only quantitative targets | 1 |
| Median validation-only factor error | 1.9340 |
| Validation-only targets within factor 2 | 100% |
| All quantitative calibration + validation targets | 29 |
| Median factor error across all quantitative targets | 1.1187 |
| All quantitative targets within factor 2 | 93.1% |
| All quantitative targets within factor 5 | 100% |

The direct/range metrics are in-sample calibration diagnostics. The M-005 value is an independent
validation of the chromium-diffusivity material property with respect to the current fit, because it
is withheld from calibration. The combined 29-target metric is a useful descriptive summary across
calibration and validation targets, but it must not be described as an independent-validation
metric.

The OpenPronghorn 0D replay is a separate software-verification test showing that the C++ model
reproduces the frozen correlation and material-property behavior.
