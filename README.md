# Aggregated HVAC System

Supplementary MATLAB code for the paper [*Modeling and control of aggregated
energy flexibility in multi-zone HVAC systems*](https://doi.org/10.1016/j.enbuild.2025.116799),
published in *Energy and Buildings*.

## Reproduce the published results

The repository includes the saved cases used by the paper. From MATLAB, run:

```matlab
addpath('scripts');
report = reproduce_paper;
assert(report.passed);
```

This command loads the reference datasets, validates the reported numerical
results, and exports the reproduced figures and validation artifacts to
`output/paper_reproduction/`.

For a numerical regression test without plotting, run:

```matlab
addpath('tests');
run_paper_validation;
```

The regression suite checks the reported aggregate parameters, prediction
errors, tracking errors, and tabulated results.

## Repository layout

- `data/published/`: immutable result data used for paper reproduction
- `scripts/`: active reproduction, analysis, and co-simulation entry points
- `src/`: reusable modeling, control, plotting, and validation code
- `tests/`: numerical regression and smoke tests
- `legacy/`: archived identification, baseline, plotting, and figure material
- `docs/`: technical notes

## Full co-simulation

The active entry points are `scripts/run_setpoint_control.m` and
`scripts/run_mass_flow_control.m`. Set `modelRoot` to a directory containing
`Threezone_buildings/` and `Fivezone_buildings/`:

```matlab
addpath('scripts');
modelRoot = 'C:/path/to/building-models';
methodCase = 'PUBLISHED';
run_setpoint_control;
run_mass_flow_control;
```

Full EnergyPlus co-simulation uses the corresponding three-zone and five-zone
model directories through `modelRoot`.

## Prerequisites

Published-result reproduction requires MATLAB (verified with R2020a). A fresh
optimization/co-simulation run additionally requires:

- MATLAB Optimization Toolbox
- YALMIP
- Gurobi Optimizer (the original project used an academic installation)
- MATLAB-EnergyPlus co-simulation toolbox (`mlep`)
- EnergyPlus and the matching building model files

## Published version

The source release associated with the paper is preserved in the
`paper-original` branch and the `v1.0-paper-original` tag.

## Citation

If you use this repository in your research, please cite the associated paper:

```bibtex
@article{hahm2026modeling,
  title   = {Modeling and control of aggregated energy flexibility in multi-zone HVAC systems},
  author  = {Hahm, Seungjun and Kim, Young-Jin},
  journal = {Energy and Buildings},
  volume  = {353},
  pages   = {116799},
  year    = {2026},
  doi     = {10.1016/j.enbuild.2025.116799}
}
```

## Contact

seungjun.hahm@postech.ac.kr

## License

See `LICENSE`.
