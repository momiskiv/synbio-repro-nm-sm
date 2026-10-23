# Focused MM32 P3 concentration pilot

Version 0.10.0 adds progress reporting and the next data-collection step. It
does not change biological rules or fit the model to a desired result.

## Question and design

Does changing the initial P3 count change the direction, timing or size of the
light response relative to no polymer? The full affinity comparison used only
1000 entities. It cannot establish the concentration dependence or whether any
loading enhances light above control.

| Condition | Initial polymers | Seeds | Duration | Work with full affinity reuse |
| --- | ---: | --- | ---: | --- |
| `no_polymer` | 0 | 2026, 2027 | 30000 a.u. | Reuse two runs |
| `p3_1000` | 1000 | 2026, 2027 | 30000 a.u. | Reuse two runs |
| `p3_3000` | 3000 | 2026, 2027 | 30000 a.u. | Two new runs |
| `p3_7000` | 7000 | 2026, 2027 | 30000 a.u. | Two new runs |

Eight trajectories are planned, with four new simulations if all four old runs
pass validation. Standard PS cutoff ten, P3 PB/PS kinetics, lattice, initial
bacteria, growth, transport and response settings are retained. Sampling remains
every 100 a.u. with 61 frames. Keeping presentation settings also permits strict
reuse. Only polymer count changes within each seed. Initial bacterial and signal
placement streams match; event histories and realised geometry can diverge.

This selected four-loading, two-seed pilot is exploratory. It is not the paper's
full sweep or ten-run average. The controls use the same two seeds as the new
conditions; ten control runs are not mixed into a two-seed paired comparison.

## 1. Activate and preview

Extract v0.10.0 to a new folder, enter `polymer-clustering`, and activate your
existing environment. Dependencies have not changed. These commands use the
extracted source directly; an editable installation is unnecessary.

Set this variable to the original completed **50-run affinity batch folder**,
containing `batch.json` and the individual condition/run folders. Replace the
placeholder. The intact original outputs on your computer are needed; an uploaded
summary JSON or the shorter crowding batch cannot supply these trajectories.

```bash
AFFINITY_BATCH="/absolute/path/to/completed/affinity/batch"

python -m polymer_clustering.batch \
  --config examples/mm32_concentration_pilot.json \
  --reuse-batch "$AFFINITY_BATCH" \
  --workers 2 \
  --dry-run
```

Expected: **8 runs; 4 reused/retained; 4 simulations pending**. Dry-run checks
completion, settings, numerical source fingerprints, NumPy, observation grids
and balances. It creates no output folders and launches no simulations.
Resolve mismatches before launching. Without reuse all eight runs are scheduled.
Changing fonts or saved-frame counts can prevent matching. Never edit original
metadata to force reuse.

## 2. Launch the four new runs

Use the same terminal so `AFFINITY_BATCH` remains defined:

```bash
python -m polymer_clustering.batch \
  --config examples/mm32_concentration_pilot.json \
  --reuse-batch "$AFFINITY_BATCH" \
  --workers 2 \
  --progress on \
  --output results/exercise5/mm32_concentration_pilot
```

Two workers are a conservative starting point while writing. More workers may
reduce responsiveness and increase memory use. Runtime is not predictable from
polymer count alone: event rates and attachment/geometry costs both change.

The batch bar starts at 4/8 after copying the four validated reused runs. It
counts verified complete jobs, including their output generation. Active jobs
show simulation percentage or their latest phase. Narrow terminals truncate the
rightmost details; optional per-attempt `progress.json` files retain each child's
latest status. Individual `console.log` files remain readable.

**Simulation 100% precedes saving plots and GIFs.** The controller counts a job
complete only after successful process exit and numerical validation. Ensemble
tables and figures are then generated, with a separate announcement. The bar
does not predict remaining real time.

Default `--progress auto` redraws in a terminal. `--progress on` also produces
throttled plain lines when redirected, without carriage returns or ANSI codes.
`--progress off` hides the live display. These options do not change model
settings, seeds or scientific results. Final runtime reporting remains.

## 3. If interrupted

Use the unique batch folder printed at launch and the same package/environment:

```bash
python -m polymer_clustering.batch \
  --resume "/actual/printed/batch/folder" \
  --workers 2 \
  --progress on
```

Completed jobs are retained; unfinished jobs replay from their seeds in new
attempt folders. Progress files are status reports, not restart checkpoints.
Event-limited jobs remain incomplete and no partial ensemble is averaged.
Cross-version import of complete runs is separate from same-version resume.

## 4. Review files and interpretation

Send these files from the printed completed batch folder:

- `batch.json`
- `ensemble.png`
- `ensemble_trajectory.csv`
- `replicate_summary.csv`
- `paired_differences.csv`
- `pairwise_differences.csv`

Keep individual logs, NPZ archives and animations locally; we can select them
later if the summary raises a spatial question. Compare mean light per cell,
then total light, receptor sensing, population and cluster size. Total light
alone can also change with population. Inspect both seeds; bands show sample
SD across two runs, not confidence intervals.

`paired_differences.csv` gives each polymer condition minus no polymer within
seed, then the mean across the two differences. Negative light differences mean
suppression and positive differences mean enhancement at that sampled time.
`pairwise_differences.csv` also compares loadings: later-listed minus earlier.

Review 5000, 15000 and 30000 a.u., sampled peak time/light and the complete curve.
Condition mean light divided by control mean light is a descriptive fold-change
when the denominator is positive. It differs from averaging per-seed ratios.
Raw light remains unchanged; publication normalisation is not imposed.

## Source audit and next decision

The supplied supplement's S26 modelled-luminescence section explicitly states
initial polymer values 0-7000 in steps of 1000 and time 0-30000 every 100.
Some legends/captions instead give 70000 or time 0-3000. This pilot follows the
explicit S26 range and selects 0, 1000, 3000 and 7000. SI S31-S39 contain MM32
time courses and fold-change surfaces; different regimes can suppress or
enhance response. Neither outcome is imposed here.

Polymer entities are a reconstruction input, without an established mg/ml
conversion. Individual on/off kinetics, initial geometry and published light
normalisation remain assumptions. Qualitative trend comparison is possible;
calibrated magnitude agreement is not established. See CONC-01 in
[model notes](model-notes.md).

Review this pilot before choosing further collection. Informative loading
effects may justify extending selected conditions to ten seeds; weak or
inconsistent effects need a targeted decision. The source's ten-run convention
is the eventual comparison target, not an automatic full sweep now.

BB170 remains a later step. Its published initial light decline needs an explicit
initial-expression diagnostic: with current `response_initial=0`, light starts
at its minimum for non-negative expression. The source does not specify E(0).
This release retains the default and does not launch a long BB170 ensemble.
Data collection and source choices will guide the eventual write-up.
