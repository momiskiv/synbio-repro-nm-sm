# Independent replicates and batch execution — v0.6.0

The next Exercise 5 question is whether the MM32 polymer/control difference is
consistent across initial conditions and stochastic histories. One trajectory
per condition cannot answer that question. This release adds batch execution and
transparent summaries; it does not change the biological reconstruction.

## 1. Install and check the plan

Run these commands from the extracted `polymer-clustering` directory, using your
existing environment or a newly activated one. The dependencies are unchanged:

```bash
python -m pip install -r requirements-validated.txt
python -m polymer_clustering --version
python -m polymer_clustering.batch --config examples/mm32_replicates.json --dry-run
```

The version should read `0.6.0`. The dry run validates settings and fonts and lists
20 jobs: two conditions, each with ten seeds, 2026–2035. It starts no simulations
and creates no output directories. Font discovery can use Matplotlib's own cache.

The scientific configuration retains your long-run controls: 25 initial MM32
cells, 10,000 signals, a 50×50 lattice, growth, 30,000 a.u., observations every
100 a.u., 61 spatial frames, and an event cap of 100 million. The conditions
differ only in polymer count: zero versus 1,000 P3 entities. Individual P3
kinetics, the concentration proxy, transport geometry, expression, and division
interpretations remain documented in `model-notes.md`.

## 2. Run the short software check first

```bash
python -m polymer_clustering.batch --config examples/replicate_smoke.json --workers 1 --output results/exercise5/batch_check
```

This runs four small jobs and writes a combined plot. **It deliberately changes
strain, lattice, division timing, and production rate** to test the machinery
quickly. Its numbers are not MM32 results or a source-matched biological test.
Every per-run figure retains the existing changed-division-timing label. The
batch plot says “raw model outputs”; always identify this diagnostic by its
configuration and output folder, not by its numerical appearance.

It lasts only two model time units, so receptor occupancy can stay at zero and
the light curves can overlap at the BB170 baseline. That is expected for this
small scheduling/output check. It is not evidence about polymer effects.

## 3. Start the ten-replicate comparison

The conservative default is one simulation at a time:

```bash
python -m polymer_clustering.batch --config examples/mm32_replicates.json --output results/exercise5/mm32_replicates
```

If your CPU allocation and available memory allow two simulations at once:

```bash
python -m polymer_clustering.batch --config examples/mm32_replicates.json --workers 2 --output results/exercise5/mm32_replicates
```

Choose one command. Starting both creates separate, duplicate batches. The
controller prints the unique batch folder immediately, then progress as jobs
start and finish. Each job has its own `console.log`; it retains the original
single-run directory and all 32 ordinary outputs, including GIFs. The larger
number of runs also multiplies event-log storage and animation memory.

On the student's machine, the seed-2026 control took about 35 minutes and P3
about 68 minutes. Twenty runs of similar cost would total roughly 17 hours of
serial work. This extrapolation is not a schedule or a measured parallel
speed-up. Two concurrent jobs can reduce elapsed batch time when resources are
available; they do not make an individual trajectory faster. Benchmark one and
two workers with the short diagnostic before choosing a long-run worker count.

Workers cannot exceed the available process-affinity/CPU count or, when set,
`SLURM_CPUS_PER_TASK`. Memory availability is not detected automatically. The
default remains one; do not request all advertised cores without checking RAM
and your allocation. Each child has common BLAS/OpenMP thread limits set to one
to avoid multiplying native threads. These limits are recorded in `batch.json`.

### Reuse the two completed seed-2026 runs

Use the original **complete run folders** containing `run.json` and the normally
named trajectory files. Do not point to a PNG, a tar archive, a renamed attachment,
or the parent containing many run folders. Extract your archive first if needed.
Replace both example paths below with your actual original run-folder paths:

```bash
python -m polymer_clustering.batch --config examples/mm32_replicates.json --workers 2 --reuse-run /path/to/no_polymer/full_RUN --reuse-run /path/to/p3/full_RUN --output results/exercise5/mm32_replicates
```

Add `--dry-run` first to verify the matches without starting work. Matching jobs
are listed as `reuse`. The controller checks complete status, settings, all six
numerical source hashes, NumPy version, observation times and signal/polymer
balances. All settings must match except the event cap, which may differ for a
run that already completed. Older v0.5.1 run metadata is retained. The controller
copies imported folders into the batch, preserving the originals, and then
launches only the 18 outstanding runs. It never counts seed 2026 twice.

### Resume after interruption

Ctrl-C terminates active child processes and preserves completed runs. Copy the
actual printed batch-folder path into:

```bash
python -m polymer_clustering.batch --resume /path/to/printed/batch_folder --workers 2
```

The saved plan and package/environment must match. A batch accepts only one
active controller. The batch controller uses POSIX process locking and has been
validated on Linux; Windows execution is not supported by this controller.
Completed numerical records
are checked again. Unfinished, failed and interrupted jobs receive a new attempt
folder and replay from their original seed. This resumes the **batch**, not the
internal stochastic state of a half-finished trajectory. Earlier attempts remain
available for inspection.

Event-limited runs are marked incomplete and are never averaged. Repeating the
same cap may reproduce that stop. To increase it, make a copy of the batch JSON,
increase `base_parameters.max_events`, and start a new batch; import the already
complete runs with repeated `--reuse-run` flags. Do not edit `batch.json` to hide
an incomplete status or change a saved plan.

## 4. Read the ensemble outputs

Only when every planned job completes and passes the boundary checks does the
controller write an ensemble average. There is no silent omission, interpolation,
or switch to a changing number of replicates at late times.

| File | Meaning |
| --- | --- |
| `batch.json` | Full job requests, seeds, source/software fingerprints, statuses, attempts, font provenance and elapsed controller time |
| `replicate_summary.csv` | One row per completed seed/condition, final state, highest sampled mean light and its time, event count and original run timing |
| `ensemble_trajectory.csv` | Condition/time/metric with replicate count, mean, sample SD and SEM |
| `paired_differences.csv` | Same statistics for condition minus reference, calculated within each shared seed first; the first listed condition is the reference |
| `ensemble.png`, `ensemble.svg` | Mean per-cell light, total light, receptor occupancy and largest connected cluster; shaded bands show pointwise ±1 SD across seeds |

Each replicate contributes one equally weighted value at each time. Cell light
is transformed before averaging within a run, as before. Run means are then
averaged across replicates; cells from large populations are not pooled across
runs. Sample SD uses `ddof=1`; SEM is SD divided by the square root of the number
of runs. For a single seed, SD and SEM are undefined and written as blank fields.
The shading is neither across-cell variation nor a confidence interval. A lower
SD band can be negative even for non-negative quantities; it is a descriptive
symmetric band, not a feasible-value prediction interval.

Within each condition, ten different seeds produce independent realisations.
The same seeds are reused across conditions to pair initial conditions and
random streams. The two conditions are therefore **paired**, rather than twenty
mutually independent trajectories. A negative P3-minus-control light difference
means lower light with P3. Summarising the differences within seeds respects that
pairing. Adjacent time points belong to the same trajectories and are correlated;
the files do not perform hundreds of independent significance tests.

Arial is requested for individual and ensemble figures. Its actual selection or
announced fallback is recorded. A reused run preserves its original font record.
To require Arial strictly for new runs, set `strict_font` to true in the batch
base parameters, supplying `font_path` if needed. That presentation change will
prevent reuse of a run saved with different requested font settings.

Controller wall time covers the current invocation, including validation,
copying reused runs, scheduling, child work and aggregation. It is recorded separately from
per-run CPU/wall measurements. Summed original run CPU seconds include imported
work and must not be presented as CPU consumed during the current invocation.

## 5. What we will review together

First check all 20 jobs completed with ten seeds per condition. Then compare the
light trajectories and their variation: does P3 suppression recur, how variable
is its timing, and how do receptor occupancy and connected-cluster size develop?
Review both mean and total light alongside population records. The existing
seed-2026 result motivates these questions; it does not prescribe their answers.

This implements the replicate machinery for the first MM32 comparison. Affinity
conditions, BB170, sensitivity analyses, and comparisons against published
experimental/model curves remain subsequent tasks. Averaging improves the
assessment of stochastic variability; it does not resolve source ambiguities or
calibrate polymer entity counts to a physical concentration.
