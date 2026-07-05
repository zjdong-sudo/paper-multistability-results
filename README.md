# Detecting the Smallest Zero-One Reaction Networks with Multistability — Code and Data

This repository contains the code and data for the paper:

> **"Detecting the Smallest Zero-One Reaction Networks with Multistability"**
> Xiaoxian Tang, Jiandong Zhang
> School of Mathematical Sciences, Beihang University

It provides all intermediate data files and source code needed to reproduce the paper's computational results: a systematic enumeration and classification of `(3,6,3)` quadratic zero-one reaction networks that admit multistability.

## Repository structure

```
├── README.md
├── code/                          # All source code
│   ├── algorithm-1.nb             # Mathematica: Algorithm 1 (Steps 1–4) — enumeration
│   ├── generate_reff-2.c          # C: Step 5a — RREF grouping
│   ├── algorithm-3.c              # C: Step 5b — species-permutation dedup
│   ├── algorithm4.nb              # Mathematica: Steps 6–7 — Gröbner basis filtering
│   ├── algorithm4_cluster.nb      # Mathematica: Step 7 — steady-state ideal clustering
│   ├── check_injectivity.nb       # Mathematica: Step 8 — injectivity screening
│   ├── check_3_steadystates.mw    # Maple: Step 9 — RealRootClassification for 3–8 roots
│   ├── maple_log.py               # Python: Step 10 — parse Maple logs
│   ├── check_multistability.mw    # Maple: Step 11 — stability verification (fixed; see §Caveats)
│   ├── check_multistability_detail.mw  # Maple: single-network detailed stability check
│   └── extract_375_networks.py    # Python: export final (Y,N) matrix pairs
├── data/                          # Intermediate pipeline files
│   ├── 363_2-algorithm1-4.txt.gz  # Master enumeration table (530,480 rows, compressed)
│   ├── 363_2-new_binary_reff_new.txt          # RREF-group assignments (104,086 groups)
│   ├── 363_2-new-binary-indices-forgbcheck.txt # Dedup representatives (409,232 indices)
│   ├── final_equivalence_representatives.txt   # Final equivalence-class reps (349,597 indices)
│   ├── potential_multistable_indices.txt       # Injectivity survivors (245,468 indices)
│   ├── odesystem.csv.gz           # ODE systems for injectivity survivors (compressed)
│   ├── candidate_indices_1.csv    # RealRootClassification candidates — batch 1 (4 indices)
│   ├── candidate_indices_2.csv    # batch 2 (64 indices)
│   ├── candidate_indices_3.csv    # batch 3 (260 indices)
│   ├── candidate_indices_4.csv    # batch 4 (47 indices)
│   ├── check_multistability_1.txt # Stability verification output — batch 1
│   ├── check_multistability_2.txt # batch 2
│   ├── check_multistability_3.txt # batch 3
│   └── check_multistability_4.txt # batch 4
└── results/                       # Final outputs
    ├── 375_multistable_networks.json  # Machine-readable: all 375 networks as (Y,N) matrices
    └── 375_multistable_networks.txt   # Human-readable: same data with reaction notation
```

## Data file conventions

The pipeline uses a consistent indexing scheme:

- **`master_index`** — 1‑based line number in `363_2-algorithm1-4.txt`. Each row is 36 comma‑separated integers ∈ {−1,0,1}. The first 18 integers are the **reactant (Y) matrix** (3 species × 6 reactions, row‑major: `[s=1,r=1..6], [s=2,r=1..6], [s=3,r=1..6]`). The next 18 are the **stoichiometric (N) matrix** in the same order. The mass‑action ODE is `dx/dt = N · v(k,x)` where `vⱼ = kⱼ · x₁^(Y₁ⱼ) · x₂^(Y₂ⱼ) · x₃^(Y₃ⱼ)`.

- **`m` (candidate index)** — 1‑based line number in `potential_multistable_indices.txt` and `odesystem.csv`. `potential_multistable_indices.txt` line `m` gives the `master_index` that row `m` of `odesystem.csv` corresponds to. All Maple‑stage files (`candidate_indices_*.csv`, `check_multistability_*.txt`) use `m`.

- **`rank`** — 1‑based order within the final 375‑network list (ascending by `m`).

## Dependencies

- **Mathematica 14.0** (or compatible) — for enumeration, Gröbner basis, and injectivity steps.
- **Maple 2024** (or compatible) with the `RegularChains`, `ParametricSystemTools`, and `SemiAlgebraicSetTools` packages — for RealRootClassification and stability verification.
- **GCC** (or any C99 compiler) — to compile `generate_reff-2.c` and `algorithm-3.c`.
- **Python 3.9+** with `sympy` (`pip install sympy`) — for log parsing and final export. The `openpyxl` package is only needed if you want to re‑run the Maple‑log parsing step (`maple_log.py`).

## Reproduction paths

### Path A — Full reproduction (from enumeration to final 373)

Follow the paper's Table 1 and Table 2. Approximate run times on a workstation (32 cores, 64 GB RAM):

| Step | Code | Tool | Input → Output | ~ Time |
|---|---|---|---|---|
| 1–4 | `algorithm-1.nb` | Mathematica | enumeration → `363_2-algorithm1-4.txt` | ~1.5 h |
| 5a | `generate_reff-2.c` | C | `363_2-algorithm1-4.txt` → `363_2-new_binary_reff_new.txt` | ~10 min |
| 5b | `algorithm-3.c` | C | above + master table → `363_2-new-binary-indices-forgbcheck.txt` | ~20 min |
| 6 | `algorithm4.nb` | Mathematica | indices + master table → Gröbner filtering | ~1 h |
| 7 | `algorithm4_cluster.nb` | Mathematica | step6 batches → `final_equivalence_representatives.txt` | ~30 min |
| 8 | `check_injectivity.nb` | Mathematica | reps + master → `potential_multistable_indices.txt` + `odesystem.csv` | ~1 min |
| 9 | `check_3_steadystates.mw` | Maple | `odesystem.csv` → RRC log | ~4 h |
| 10 | `maple_log.py` | Python | RRC log → `candidate_indices_*.csv` | seconds |
| 11 | `check_multistability.mw` | Maple | candidates + `odesystem.csv` → `check_multistability_*.txt` | ~1 min |
| — | `extract_375_networks.py` | Python | all above → `results/*` | < 1 s |

**Before starting:** decompress the two compressed files:
```bash
gunzip data/363_2-algorithm1-4.txt.gz
gunzip data/odesystem.csv.gz
```

The C programs read input files from the current working directory and write output there as well. The Mathematica, Maple, and Python scripts use hard‑coded Windows‑style paths (`D:/multistability/…`) in the versions archived here — you will need to change these to your local paths before running.

### Path B — Quick start (from injectivity survivors onward)

If you only want to re‑run the Maple multistationarity / multistability detection (Steps 9–11), you can skip the expensive Gröbner‑basis and injectivity steps and start directly from the provided `odesystem.csv`:

1. Decompress: `gunzip data/odesystem.csv.gz`
2. Run `check_3_steadystates.mw` (Maple) on `odesystem.csv` — this runs RealRootClassification on all 245,468 injectivity‑surviving networks with a 60‑second timeout per network.
3. Run `maple_log.py` (Python) to parse the Maple log and extract candidate indices.
4. Run `check_multistability.mw` (Maple) on the candidates to verify stability with numerical sample points.
5. Run `extract_375_networks.py` (Python) to produce the final `(Y,N)` matrix pairs and an audit summary.

## Final results

The file `results/375_multistable_networks.json` contains all 375 RRC‑surviving networks as `(Y,N)` matrix pairs with derived reaction equations, ODE self‑consistency checks, and per‑sample‑point stability data. The `_audit_summary` key at the top of the JSON reports:

- **373 networks** confirmed multistable (≥ 2 stable positive steady states in at least one parameter sample point).
- **2 networks flagged** (`m=33449`, `m=239089`): the RRC step correctly found 3 positive steady states for each, but at every numerically‑solved sample point only 1 of the 3 is stable. See the `multistability_confirmed` and `per_sample_summary` fields on those two records for details.

All 544 sample points across all 375 networks produced exactly **3** positive steady states — none showed 4 or more, despite the Bézout upper bound of 8.

## Caveats

- **Path dependency.** The Mathematica and Maple scripts contain hard‑coded Windows file paths (`D:/multistability/…`). Adjust these to your local directory structure before running.

- **Non‑determinism in parallel steps.** The Mathematica `ParallelTable` calls and Maple's sample‑point selection can produce slightly different outputs (e.g., a different numeric sample point or a different ordering of Gröbner‑basis results) across runs and platforms. The equivalence‑class deduplication and the final count should be invariant; individual row indices may differ.

- **The 239089 and 33449 networks.** The original Maple stability‑verification step (before the fix described below) accumulated stable‑steady‑state counts across *different* sample points in its summary table, which could make a network appear to have `stable ≥ 2` when no single sample point individually reached that threshold (see `m=33449`: two sample points, each with 1 stable, summed to `stable=2`). The current version of `check_multistability.mw` in this repository has been corrected to report the per‑sample‑point maximum and a boolean `is_multistable` flag. The extraction script `extract_375_networks.py` independently computes `multistability_confirmed` from per‑sample‑point detail blocks, bypassing the summary table entirely.

- **RealRootClassification range.** RRC was called with range `3..8`. It asks "does there exist a parameter region where the system has *some* number in [3,8] of positive real roots?" — it does not distinguish between "exactly 3" and "4, …, 8". The empirical result that all 375 networks admit *exactly* 3 positive real roots is based on numerically solving the equations at all 544 sample points; no symbolic proof of this exact‑3 upper bound is claimed.

## License

This repository is provided to accompany the paper. Please cite the paper if you use this code or data in your own work.
