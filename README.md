# EMILI IMRT — Beam Angle & Fluence Map Optimization

<p align="center">
  <pre>
   ______ __  __ _____ _      _____ 
  |  ____|  \/  |_   _| |    |_   _|
  | |__  | \  / | | | | |      | |  
  |  __| | |\/| | | | | |      | |  
  | |____| |  | |_| |_| |____ _| |_ 
  |______|_|  |_|_____|______|_____|
  </pre>
</p>

**EMILI** (*Easily Modifiable Iterated Local search Implementation*) is an extensible C++ metaheuristic framework for combinatorial optimization, extended here with an **IMRT treatment planning module** that solves the coupled **Beam Angle Optimization (BAO)** and **Fluence Map Optimization (FMO)** problem using AMPL + Gurobi and iterated local search.

---

## Table of Contents

- [Architecture](#architecture)
- [IMRT Problem Formulation](#imrt-problem-formulation)
- [Quick Start](#quick-start)
- [Command-Line Usage](#command-line-usage)
- [Hyperparameter Tuning (irace)](#hyperparameter-tuning-irace)
- [Clinical Report (ICRU-83)](#clinical-report-icru-83)
- [Algorithm Competition](#algorithm-competition)
- [Project Structure](#project-structure)
- [Build Options](#build-options)
- [Contributing](#contributing)
- [License & Credits](#license--credits)

---

## Architecture

The codebase follows a **framework + extension** pattern. The core EMILI library provides abstract metaheuristic building blocks; the `imrt/` directory extends them for radiation therapy.

```
┌──────────────────────────────────────────────────┐
│  FRAMEWORK CORE  (emilibase.h/.cpp)              │
│  Problem · Solution · InitialSolution            │
│  Neighborhood · Perturbation · Acceptance        │
│  Termination · LocalSearch · TabuMemory          │
├──────────────────────────────────────────────────┤
│  PARSER / BUILDER  (generalParser.h/.cpp)        │
│  CLI grammar → algorithm composition at runtime  │
├──────────────────────────────────────────────────┤
│  IMRT MODULE  (imrt/)                            │
│  CerrFmoSource · ImrtProblem · ImrtFmoSolver     │
│  BaoProblem · AngleSwapNeighborhood              │
│  RandomAnglesPerturbation · GreedyAngles...      │
└──────────────────────────────────────────────────┘
```

### Metaheuristics available out of the box

| Algorithm | Token | Internal LS | Special components |
|---|---|---|---|
| Iterated Local Search | `ils` | Best/First improvement | `prangswap`, `improve` |
| Variable Neighborhood Search | `vns` | Best/First improvement | `bangshake`, `accng` |
| Tabu Search (BAO) | `tabubao` | none (own loop, neighbours from the current solution) | `TBao_fixed`, `TBao_adaptive`, `TBao_attr` |
| Variable Neighborhood Descent | `vnd` | list of neighborhoods | `nangshift`, `nangswapsample`, `nangswapwindow`, `nangflip` |
| Iterated Greedy | `ils` + `pgreedy` / `pgreedysample` | Best/First improvement | destruction-reconstruction (sampled variant) |
| Simulated Annealing | `ils nols` + `saacc` | Empty (Metropolis accept) | `saacc s e β i α` |

*All algorithms are composed at runtime via CLI tokens — no recompilation needed.*

---

## IMRT Problem Formulation

### Clinical context

Intensity-Modulated Radiation Therapy delivers conformal dose to a Planning Target Volume (PTV) while sparing Organs at Risk (OARs). The problem has two coupled sub-problems:

1. **BAO (Beam Angle Optimization)** — select *K* gantry angles from *N* candidates
2. **FMO (Fluence Map Optimization)** — compute optimal beamlet intensities for the chosen angles

### Quadratic penalty objective (FMO)

The FMO sub-problem is solved **exactly** via AMPL + Gurobi:

```
min  w_under · Σ || u_b ||²  +  w_over · Σ || v_b ||²  +  w_ptv_over · Σ || w_b ||²
s.t. D_ptv·x + u ≥ Dmin            (PTV underdose slack)
     D_oar·x  − v ≤ Dmax           (OAR overdose slack)
     D_ptv·x  − w ≤ 1.07·Dmin      (PTV hot-spot constraint)
     u, v, w ≥ 0
     0 ≤ x_j ≤ max_intensity
```

| Weight | Controls | Clinical effect of increasing |
|---|---|---|
| `w_under` | PTV underdose penalty | ↑ forces V95% coverage |
| `w_over` | OAR overdose penalty | ↑ protects healthy organs |
| `w_ptv_over` | PTV hot-spot penalty | ↓ reduces D2 / improves homogeneity |

### Sparse dose matrix

Dose deposition data is read from the raw **CERR export layout** (`instances/CERR_Prostate`) via `CerrFmoSource` — one voxel list per organ plus per-angle dose files, with beamlet counts read per angle instead of assumed uniform. The FMO solver assembles a **compact QP** using only the *K* active angles' beamlets, keeping the problem size O(K) regardless of candidate pool size.

---

## Quick Start

### Prerequisites

- C++11 compiler (gcc ≥ 4.7 or Clang)
- CMake ≥ 2.8
- AMPL + Gurobi, installed locally via `ampl_gurobi/.venv` (see `ampl_gurobi/` for setup — the AMPL C++ API itself is vendored under `third_party/ampl_cppapi/`, no separate install needed for that part)

### Build

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The FMO solver (`imrt/imrt_fmo.cpp`) resolves the AMPL binary and Gurobi
driver at runtime via `EMILI_AMPL_BIN_DIR`/`EMILI_GUROBI_BIN` (defaulting to
the paths inside `ampl_gurobi/.venv`) — if that venv isn't set up, `baoimrt`
and `imrt` will still build and run, but every FMO solve will fail loudly
and return a worst-case objective instead of a real result.

The binary is `build/emili`.

---

## Command-Line Usage

```bash
./build/emili <instance_dir> <problem> <algorithm_description> [options]
```

### BAO+FMO (K=4 angles, ILS with first improvement)

```bash
./build/emili instances/PROSTATE_36ang baoimrt 4 \
    ils first irandomk tmaxiter 1 nangswap \
    tmaxiter 3 prangswap 2 improve \
    rnds 42
```

### FMO-only (all angles, single-beamlet-shift neighborhood)

```bash
./build/emili instances/TINY_test imrt 36 \
    bi izero nshift 0.5 tmaxiter 100 \
    rnds 42
```

### Algorithm tokens quick reference

| Category | Token | Meaning |
|---|---|---|
| Problem | `imrt <dir> [nactive K] [verbose]` | FMO on all angles |
| Problem | `baoimrt <K> <dir> [verbose]` | BAO+FMO |
| Init | `izero` / `iuniform X` / `irandom M` | FMO initial intensities |
| Init | `ifirstk` / `irandomk` | BAO initial angle set |
| Neighborhood | `nshift Δ` / `nswap` | FMO beamlet moves |
| Neighborhood | `nangswap` | BAO angle swap |
| Perturbation | `prangswap P` / `pgreedy D` | BAO perturbations |
| Acceptance | `improve` / `saacc s e β i α` | Accept criterion |
| Termination | `tmaxiter N` / `tfeasible` | Stop condition |

---

### BAO token reference

Neighborhoods (all change one angle of the K active angles):

| Token | Neighbours | Notes |
|---|---|---|
| `nangshift <step> [fixed\|circular\|randorder]` | 2K (angle -step, +step) | step in degrees (1, 5, 10, 20); the scan order only matters for first improvement, default `fixed` |
| `nangswap` | K x (catalog - K) | full swap, expensive (272 neighbours with the 5-degree catalog) |
| `nangswapsample <n> [unrestricted]` | n sampled swaps | a plan with no improving sample is only probably a local optimum |
| `nangswapwindow <deg> [unrestricted]` | swaps within +-deg of the removed angle | |
| `nangflip 0` | at most K (angle +180) | |

Perturbations (ILS and the `pershake` shake of the general VNS):

| Token | Meaning |
|---|---|
| `prangshift <step> <n> [stepwise] [distinct\|replace]` | shifts n angles by a random magnitude in (step, 2 step); `distinct` (default) moves n different angles, `replace` draws the slot of each move independently, so an angle can move more than once; `stepwise` evaluates after every individual move |
| `prangshiftadaptive <step> <max> <stagnation> [distinct\|replace]` | starts at 1 angle and adds one per `stagnation` rounds without improvement |
| `prangshiftmag` / `prangshiftadaptivemag` | same, with discrete magnitudes |
| `pgreedy <D>` / `pgreedysample <D> <sample>` | destroy D angles and rebuild greedily (the sampled variant evaluates `sample` candidates per empty slot) |

Tabu search:

```
tabubao [first|best] [stall <n>] irandomk tmaxiter <N> nangshift <step> <memory>
```

- `first|best` is optional (default `best`, the textbook strategy); `stall` stops after n iterations without a new FMO solve (default 100).
- Memories: `TBao_fixed <tenure>` and `TBao_adaptive <min> <max>` forbid whole angle sets already visited; `TBao_attr <tenure>` forbids re-inserting any of the last `tenure` angles that left the plan.
- Debug output: `EMILI_TABU_LOG=<file>` writes one csv row per tabu event (start, move, tabu_skip) with the banned list before and after; `EMILI_PERT_TRACE=1` prints the slot and angle change of every perturbation move.

Example:

```bash
./build/emili instances/CERR_Prostate baoimrt 4 csv trajectory.csv \
    tabubao first irandomk tmaxiter 1000000 nangshift 10 TBao_adaptive 3 8 \
    -it 1080 rnds 1
```

`-it <seconds>` must come before `rnds`. The XML grammar `grammars/imrt_bao.xml` generates these command lines for automatic configuration; `scripts/check_irace_params.py` checks it against the irace parameter space and `scripts/sample_grammar.py` prints random samples.

---

## Per-organ FMO model and experiment variables

The default FMO model is `ampl_gurobi/fmo.mod` (one overdose weight for all organs at risk). The per-organ model
`ampl_gurobi/fmo_organs.mod` gives each organ its own weight. It is selected at run time, without recompiling:

| Variable | Default | Meaning |
|---|---|---|
| `EMILI_FMO_MODEL` | `baseline` | `baseline` reads `fmo.mod`; `organs` reads `fmo_organs.mod` |
| `EMILI_W_OVER_ORGANS` | `w_over` for every organ | Overdose weights, comma separated, in the order of the organs at risk (e.g. `0.5,0.8`) |
| `EMILI_FMO_NORMALIZE` | `0` | `1` divides each term by the number of voxels of its structure |
| `EMILI_DMIN` | `65` | PTV prescription in Gy (the PTV ceiling is 1.07 x this value) |
| `EMILI_W_PTV_OVER` | `0` | Weight of the PTV hot-spot penalty |
| `EMILI_EXCLUDE_OVERLAP` | `0` | `1` counts voxels shared by a PTV and an organ only as PTV |
| `EMILI_MAX_INTENSITY` | `15000` | Upper bound of each beamlet intensity |
| `EMILI_GUROBI_OPTIONS` | empty | Passed to the solver, e.g. `threads=1` |
| `EMILI_AMPL_BIN_DIR`, `EMILI_GUROBI_BIN` | repo `.venv` path | AMPL and solver binaries (needed on Linux) |

Defaults reproduce the original behaviour. The runners in `scripts/` set the per-organ variables used in the experiments
(for example `scripts/run_catalog_job.sh`). The solver license is not part of the repository (`*.lic` is ignored);
see `ampl_gurobi/README.md`.

## Hyperparameter Tuning (irace)

`irace/` holds the [irace](https://iridia.ulb.ac.be/irace/) setup used to configure the BAO metaheuristics:

| File | Content |
|---|---|
| `irace/parameters.irace` | parameter space: algorithm, strategy, step, scan order, neighbourhoods, perturbations, tabu tenures, cooling |
| `irace/scenario.irace` | budget and instances |
| `irace/target-runner` | builds the EMILI command line from the sampled parameters and returns the objective |
| `irace/results/irace.Rdata` | final result of the executed run |

```bash
irace --scenario irace/scenario.irace
```

The grammar `grammars/imrt_bao.xml` is a superset of this space (steps 1 and 20, `replace` slots); the differences are
documented in `scripts/check_irace_params.py`.

---

## Clinical Report (ICRU-83)

Every BAO+FMO run outputs a clinical plan summary:

```
================================================================
  IMRT FMO RESULT   K=4 beams
================================================================
Selected angles (deg) : [ 160, 220, 270, 290 ]
FMO objective f*      : 231180.90
Prescription dose Rx  : 68.00 Gy

----------------------------------------------------------------
  Per-organ dose statistics  (Gy)
----------------------------------------------------------------
Structure    Voxels    Dmin   Dmean    Dmax     D95      D5      D2
----------------------------------------------------------------
PTV_68         6770   54.28   67.53   79.16   60.77   73.83   74.92
Bladder       11596    0.13   35.29   70.02       -   63.10   64.21
Rectum         1764    0.04   38.55   67.51       -   62.20   64.70

----------------------------------------------------------------
  DVH constraints  (clinical goals)
----------------------------------------------------------------
Structure  Metric          Goal        Achieved      Status
----------------------------------------------------------------
PTV_68     D95 >= 0.95*Rx  64.60       60.77         VIOL
PTV_68     D2  <= 1.07*Rx  72.76       74.92         VIOL
Bladder    Dmax <= limit   50.00       70.02         VIOL

----------------------------------------------------------------
  Plan quality indices
----------------------------------------------------------------
Conformity Index    CI (V_Rx / V_PTV)    = 1.001
Homogeneity Index   HI (D2 - D98)/D_Rx   = 0.229
PTV coverage        V95% of Rx           = 71.85 %
```

A **DVH CSV** is also generated (`dvh.csv`). The search trajectory of a run (`csv <file>` option) can be plotted with `scripts/plot_trajectory.py`.

### Plan quality metrics

| Metric | Definition | Target |
|---|---|---|
| **V95%** | % PTV volume receiving ≥ 95% Rx | ≥ 95% |
| **CI** | Conformity Index: V_Rx / V_PTV | → 1.0 |
| **HI** | Homogeneity Index: (D2 − D98) / Rx | → 0 |
| **D95** | Dose to 95% of PTV | ≥ 0.95·Rx |
| **D2** | Near-max dose (2% volume) | ≤ 1.07·Rx |

---

## Algorithm Competition

`run_comparison.sh` benchmarks multiple algorithms with identical seeds:

```bash
./run_comparison.sh
```

Output per algorithm in `comparison_results/`:

| Algorithm | Log file | Key components |
|---|---|---|
| ILS | `ILS.log` | ILS + random init + 2-swap perturbation |
| Tabu | `Tabu.log` | Tabu with BAO angle-set memory |
| VNS | `VNS.log` | rVNS with multi-scale shake |
| SA | `SA.log` | ILS + Metropolis cooling |
| Greedy | `Greedy.log` | ILS + destruction-reconstruction |

Trajectory and comparison plots are in `scripts/` (`plot_trajectory.py`, `plot_aggregate_comparison.py`, `plot_boxplot_comparison.py`).

Per-seed analysis of the tabu runs: `scripts/tabu_replay_events.py` rebuilds the tabu-list events from each `trajectory.csv`
(validating the replay against the logged evaluations) and `scripts/tabu_seed_convergence_csv.py` writes `convergence.csv`
(best objective so far against real FMO solves).

---

## Project Structure

```
emili_imrt/
├── main.cpp                   # Entry point
├── emilibase.h/.cpp           # Core framework (Problem, Solution, LS, ILS, VNS, Tabu...)
├── generalParser.h/.cpp       # CLI token parser + Builder pattern
├── setup.h                    # Compile-time flags
├── CMakeLists.txt             # Build system
│
├── imrt/                      # IMRT extension module
│   ├── cerr_instance.h/.cpp   # CerrFmoSource: raw CERR export loader
│   ├── imrt_fmo_source.h/.cpp # IFmoDataSource interface shared by FMO/BAO
│   ├── imrt.h/.cpp            # FMO: ImrtProblem, beamlet neighborhoods
│   ├── imrt_fmo.h/.cpp        # AMPL+Gurobi FMO solver
│   ├── imrt_bao.h/.cpp        # BAO: angle neighborhoods, perturbations
│   └── imrt_builder.h/.cpp    # CLI builder for IMRT tokens
│
├── grammars/imrt_bao.xml      # Grammar of the BAO algorithms (grammar2code format)
├── irace/                     # irace configuration (parameters, scenario, target-runner)
├── ampl_gurobi/               # FMO models, solver drivers and dosimetry tools
├── scripts/                   # Experiment runners, analysis and plotting tools
├── experiments/               # Output of the experiment runners (not distributed)
└── run_comparison.sh          # Multi-algorithm benchmark runner
```

---

## Build Options

| CMake flag | Default | Description |
|---|---|---|
| `DEBUG_FLAGS` | ON | Debug symbols (`-g`) |
| `O3_FLAGS` | OFF | `-O3` optimization |
| `IRACE_OPTIMISED_FLAGS` | OFF | irace-tuned GCC flags |
| `OPT_OPTIMISED_FLAGS` | OFF | OpenTuner-tuned GCC flags |
| `WITH_STATS` | ON | Solution statistics output |
| `USE_NEW_MAIN` | ON | Component-based algorithm loading |
| `PORTABLE` | OFF | Static linking |

---

## Contributing

This repository follows a **PR-based workflow** managed by Gentle AI. All changes go through:

1. **Feature branch** off `develop`
2. **Code changes** with conventional commits
3. **Pull Request** with automated review
4. **Merge to `develop`** after approval

### Conventional commit prefixes

- `feat:` — new feature
- `fix:` — bug fix
- `perf:` — performance improvement
- `refactor:` — code restructuring
- `docs:` — documentation
- `ci:` — CI/irace configuration

### Instance data

Clinical instance files (Dij matrices, VOILists) are **not tracked** in this repository due to size. Contact the maintainer for access to training instances.

---

## License & Credits

- **EMILI framework**: Federico Pagnozzi (BSD 2-Clause) — `federico.pagnozzi@ulb.ac.be`
- **IMRT extension** (BAO+FMO, AMPL+Gurobi integration, irace tuning): Marco Rojas
- **AMPL / Gurobi**: AMPL Optimization Inc. / Gurobi Optimization, LLC

```
EMILI IMRT — BAO & FMO optimization for radiation therapy planning
Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
This file is distributed under the BSD 2-Clause License. See LICENSE.txt
```
