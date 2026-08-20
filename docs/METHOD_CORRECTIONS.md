# Paper-method corrections

The `paper-original` branch and `v1.0-paper-original` tag preserve the code
that generated the published results. The `PUBLISHED` experiment configuration
and regression tests keep that implementation separate from correction cases
A-D.

Figure 4 generation is intentionally outside the reproduction scope. Thermal
parameter identification remains an upstream requirement because its outputs
feed the remaining figures and Tables 5-8.

## Confirmed differences

1. Equation (30) includes
   `abs(alpha_VESS - alpha_n) / beta_VESS`, but the original implementation
   omits this correction from the aggregate charging and discharging limits.
2. Equations (13) and (38) define building SOC as a `delta_i / b_i` weighted
   average. The mass-flow controller uses an unweighted mean for online state
   estimation and for the local balancing constraint. Post-processing already
   uses the weighted definition.
3. The paper constrains SOC to `[-1, 1]`. The published implementation used
   `[-5, 5]`.

Diagnostics on the saved published cases show that these are not equivalent
implementations, even though the first-order numerical impact is modest:

- Equation (30) changes the mean discharge limit from 73.3902 to 73.1710 kW
  and the mean charge limit from 25.4520 to 24.9085 kW.
- Weighted versus unweighted mass-flow building SOC differs by 0.022908 MAE,
  0.034067 RMSE, and 0.123843 maximum absolute error.
- The published regulation signal ranges from -0.929824 to 0.935816, so the
  original `[-5, 5]` constraint is inactive in the saved case. This does not
  prove that the wider bound is harmless for other signals.

Table 5 also contains a presentation discrepancy: the published value for
VB 18 charging power is 10.6 kW, while the saved case evaluates to
10.671936 kW (10.67 kW at two decimal places). The regression test records
this 0.07 kW difference explicitly.

## Required correction experiment

Apply each correction separately before combining them:

- A: Equation (30) limit correction only.
- B: Weighted building SOC in mass-flow state estimation and optimization.
- C: SOC bounds changed from `[-5, 5]` to `[-1, 1]`.
- D: A, B, and C together.

For every experiment, compare Figures 5-13 and the metrics emitted by
`reproduce_paper`. Figure 4 remains excluded.

The diagnostic values above are formula evaluations on saved data, not results
from rerunning the corrected closed-loop co-simulation. A defensible corrected
result set requires the external EnergyPlus building models and each A-D run.

## Implementation status

- `vess_power_limits` implements Equation (30) and is shared by both main
  simulations and the saved-data diagnostics.
- `aggregate_zone_soc` and `zone_soc_weights` implement Equations (13) and
  (38) in online state estimation and mass-flow allocation.
- `schedule_demand_response` accepts the experiment-specific symmetric SOC
  bound.
- `method_correction_config` defines `PUBLISHED` and correction cases A-D.
- Simulations write to `output/method_corrections/` with names such as
  `mass_flow_control_a.mat`; published cases are not overwritten.

Set `methodCase` before running either main script. If it is omitted, case D
is used:

```matlab
addpath('scripts');
modelRoot = 'C:/path/to/building-models';
methodCase = 'A';
run_mass_flow_control;
```

`PUBLISHED` reproduces the original code choices and is not a claim that those
choices match all equations in the paper.

Run `tests/run_method_correction_tests.m` for configuration and formula-level
tests. This does not replace the A-D EnergyPlus co-simulation experiment
described above.

Once the setpoint and mass-flow outputs for a case are available, run:

```matlab
addpath('scripts');
report = analyze_method_case('A');
```

The analyzer compares 27 scalar metrics with the checked-in published cases,
writes a CSV/MAT summary, and generates Figures 5-13. It does not generate
Figure 4. `tests/run_method_analysis_smoke_test.m` verifies this pipeline by
comparing the published cases with themselves.
