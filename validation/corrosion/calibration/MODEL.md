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
and a reduced empirical magnitude. Additional closures map the common corrosion/deposition drive
to:

- uniform and intergranular attack depth;
- areal mass loss and gain;
- salt chromium and iron inventory changes;
- chromium solid-state diffusivity;
- noble-metal off-gas fraction; and
- tellurium solubility and redox threshold.

These closures are jointly fitted effective responses. Agreement in one branch does not establish a
mechanistic validation of every underlying process.

## Relation-aware residuals

The calibration does not convert every literature statement into an artificial exact value:

- Direct and range targets use weighted log residuals.
- Upper and lower bounds use one-sided residuals that are zero when the inequality is satisfied.
- The deposition-ranking row contributes two ordering residuals.
- Input-only and excluded auxiliary rows remain in the audit tables but do not enter the objective.

For a direct target,

```text
r_i = sqrt(quality_weight_i)
      * ln(prediction_i / target_i)
      / sigma_ln_i.
```

The optimizer minimizes a robust `soft_l1` loss over the data residuals and 61 Gaussian-style prior
residuals:

```text
r_prior,j = prior_strength * (parameter_j - prior_j) / prior_sigma_j.
```

The default `prior_strength` is `0.35`. Bounds, initial values, prior centers, and prior widths are
defined next to every parameter in `src/msr_corrosion_bv/model.py`.

## Optimization

The fit uses `scipy.optimize.least_squares` with:

- bounded parameters;
- `loss="soft_l1"`;
- `f_scale=1`;
- Jacobian-based variable scaling;
- a default maximum of 4000 function evaluations; and
- deterministic initial values.

The frozen run terminated successfully after 24 function evaluations. Because the fit is
deterministic, no random seed is involved.

## Interpretation of reported performance

The current fit used 43 measurement rows, 38 active constraint rows, and 61 fitted parameters plus
priors. Its frozen metrics were:

| Metric | Value |
|---|---:|
| Direct/range targets | 28 |
| Median direct/range factor error | 1.117 |
| Direct/range targets within factor 2 | 92.9% |
| Direct/range targets within factor 5 | 100% |
| Active constraint pass fraction | 92.1% |

These are in-sample calibration diagnostics. They are not an independent test-set accuracy claim.
The OpenPronghorn 0D replay is a separate software-verification test showing that the C++ model
reproduces this correlation.
