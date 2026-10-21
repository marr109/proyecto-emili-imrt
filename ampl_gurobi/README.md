# ampl_gurobi — which model is which

| File | Role |
|---|---|
| `fmo_organs.mod` | **Official model for the Python/AMPL test pipeline.** Per-organ OAR overdose weights (`w_over_organ`) and optional `normalize`. Every script here (`cerr_fmo_ampl.py`, `fi_search_*.py`, `pareto_sweep.py`, `dosimetry_report.py`) loads it. |
| `fmo_organs_envio.mod` | Short, comment-trimmed copy of `fmo_organs.mod` (identical code) for sharing. |
| `fmo.mod` | **Baseline (single `w_over`).** Frozen, do not edit: it keeps earlier experiments reproducible. |

`fmo.mod` stays in this folder on purpose: `imrt/imrt_fmo.cpp:137` loads
`ampl_gurobi/fmo.mod`, so moving it would break the C++ solver. With equal
`w_over_organ` and `normalize = 0`, `fmo_organs.mod` returns the same objective
as `fmo.mod`.

Not wired yet: `ImrtFmoSolver` (C++) does not pass `oar_organ` / `w_over_organ`.

Loader options (environment variables, `cerr_fmo_ampl.py`):
`EMILI_DMIN`, `EMILI_W_PTV_OVER`, `EMILI_EXCLUDE_OVERLAP=1` (OAR minus PTV).

## Local (macOS) Python environment

Persistent venv outside the project and outside /tmp (macOS cleans old files in /tmp and broke the previous one):
`~/.venvs/emili-ampl` (numpy, scipy, matplotlib, amplpy, ampl_module_base, ampl_module_gurobi).

```bash
~/.venvs/emili-ampl/bin/python dosimetry_batch.py <experiments_dir> --out <csv> --geud
~/.venvs/emili-ampl/bin/python ../scripts/analyze_new_experiments.py
```

The AMPL course license uses short-term leases. If a run fails with "Short term lease expired", renew it:

```bash
cd ~/.venvs/emili-ampl/lib/python3.14/site-packages/ampl_module_base/bin && ./amplkey renew
```
