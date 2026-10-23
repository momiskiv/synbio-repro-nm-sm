# Polymer clustering: Exercise 5 reconstruction

An educational Python reconstruction of the spatial model in **Lui et al.
(2013)**, DOI [10.1038/nchem.1793](https://doi.org/10.1038/nchem.1793).

**Version 0.11.1 adds the BB170 collection plans for Exercise 5 Part 3.**
Run them on the student's laptop using the optimised solver. Start with a
full-duration no-polymer/P3 pair, then reuse it in the standard eight-loading,
ten-seed study. Separate 80-run plans cover lower density and faster growth.
The faster-growth condition now accepts a zero lower endpoint in the
continuous triangular lifetime distribution while keeping sampled lifetimes
positive. Existing positive-support trajectories preserve their random draws.

From the extracted project with the validated environment active:

```bash
python -m polymer_clustering --version
python -m polymer_clustering.batch --config examples/bb170_standard_pair.json --workers 2 --dry-run
python -m polymer_clustering.batch --config examples/bb170_standard_pair.json --workers 2 --progress on --output results/exercise5/bb170_standard_pair
```

Expect version 0.11.1 and two full 30000-a.u. trajectories, seed 2026.
The [BB170 walkthrough](docs/bb170-walkthrough.md) gives the complete 240-run
source-based plan, reuse commands, interpretation choices and review files.
No new dependencies are needed. Paper comparisons follow data collection.

**Version 0.11.0 removes an expensive global cache rebuild after each polymer
bond change.** Targeted dependency updates preserve the old seeded event/time
arithmetic. Three-repeat bounded tests measured 6.9x/14.1x/25.3x solver speedups
at 1000/3000/7000 polymers. These are short diagnostic measurements, not full-run
runtime predictions. Completed scientific outputs remain reusable after the
existing validation. No dependency update is needed.

Full mode couples binding, movement, expression, predicted light, division and
strain-specific signal production. Existing affinity/crowding controls, spatial
analysis, stable-scale animations, Arial-first plots, timings, parallel batches,
resume/reuse controls and diagnostic modes remain available.

This is an implementation milestone. Source ambiguities remain documented, and
agreement with the publication's time courses has not yet been established.

## Review the optimisation (v0.11.0)

See the [performance review](docs/performance-review.md) for the profile,
implementation, exact numerical checks, raw benchmark evidence, reproduction
commands and remaining limits. Extract this release to a new folder and use it
for future runs; retain original packages for resuming their unfinished batches.
Your completed concentration pilot does not need to be repeated.

A short software check from the new project folder is:

```bash
python benchmarks/profile_solver.py --case trajectory --polymers 7000 --events 5000 --repeats 3 --profile
```

The 5000-event cap is intentional software verification. It is not a completed
scientific trajectory. Keep the old package to compare the same command using
`--source-root "/absolute/path/to/v0.10.0/polymer-clustering"`.

## Reproduce the concentration pilot (introduced in v0.10.0)

The commands below reproduce the collection plan from the earlier affinity
batch. Existing complete concentration outputs can instead supply all eight
jobs via `--reuse-batch`; a dry-run then expects eight reused and zero pending.

Set `AFFINITY_BATCH` to your original completed 50-run affinity batch folder,
containing `batch.json` and the individual run folders. Preview first:

```bash
AFFINITY_BATCH="/absolute/path/to/completed/affinity/batch"
python -m polymer_clustering.batch --config examples/mm32_concentration_pilot.json --reuse-batch "$AFFINITY_BATCH" --workers 2 --dry-run
```

Expect **4 reused/retained and 4 pending** out of eight jobs. Then launch:

```bash
python -m polymer_clustering.batch --config examples/mm32_concentration_pilot.json --reuse-batch "$AFFINITY_BATCH" --workers 2 --progress on --output results/exercise5/mm32_concentration_pilot
```

Conditions are no polymer and P3 at 1000/3000/7000 entities, seeds 2026/2027,
to 30000 a.u. with standard PS cutoff ten. This is an exploratory pilot. See the
[concentration walkthrough](docs/concentration-walkthrough.md) for reuse,
progress, interruption, review files, source choices and interpretation.
Default `--progress auto` redraws in a terminal; `on` also writes plain periodic
lines when redirected; `off` hides progress. Single-run commands accept the same
option. Progress options do not belong in scientific JSON settings.

## Test PS crowding with bacterial binding retained (v0.9.0)

The pilot compares standard P3 with P3 whose PS-site cutoff is disabled, using
two matched seeds to 10000 a.u. PB/PS kinetics, polymer count, transport and
source-default growth settings remain identical. Preview and launch four jobs:

```bash
python -m polymer_clustering.batch --config examples/mm32_crowding_pilot.json --workers 2 --dry-run
python -m polymer_clustering.batch --config examples/mm32_crowding_pilot.json --workers 2 --output results/exercise5/mm32_crowding_pilot
```

This is a new simulation pilot. See the [crowding walkthrough](docs/crowding-walkthrough.md)
for the design, interpretation, review files, limits and legacy compatibility.

## Analyse saved polymer crowding (v0.8.0)

Use completed affinity results directly, without rerunning simulations:

```bash
python -m polymer_clustering.spatial \
  --batch "/path/to/completed/affinity/batch" \
  --conditions p3_no_pb_1000 p3_1000 \
  --map-seed 2026 \
  --output results/exercise5/spatial_affinity
```

Replace the placeholder with the folder containing the completed `batch.json`.
The command writes to a separate unique folder and prints elapsed time and the
`spatial_summary.zip` to upload. Outputs include per-seed and paired endpoint
metrics, shared-scale polymer/crowding/free-signal maps, and provenance. Existing
frames lack polymer coordinates over time, so this reports endpoints rather
than inferring a crowding trajectory or causal mechanism. See the
[spatial walkthrough](docs/spatial-walkthrough.md) for preview commands,
metric definitions, scientific limits and output interpretation.

## Install and run on Linux

Requires Python 3.11 or later. From this project's directory, use your existing
virtual environment, or create one:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-validated.txt
python -m polymer_clustering
```

The default run needs no settings edits. It places 25 bacteria on a 50 x 50 grid
with 10000 initially free signals (MM32), no polymers, and seed 2026. The paper
does not prescribe our initial placement algorithm or seed; see model notes.

If your environment cannot install these pinned versions, the supported ranges
are in `pyproject.toml`. A fresh dependency resolution has not been validated here.

## Execution time and plot fonts (v0.5.1)

Every finished run prints real elapsed seconds for the total run, numerical
simulation, and output generation. Detailed timings, CPU seconds, and solver
throughput are saved under `timing` in `run.json`, including event-limited runs
and caught runtime failures. For an initial-state-only run, simulation time is
zero. These real seconds are separate from the model's arbitrary time units.

For example, the terminal format is `Runtime: total ... s; simulation ... s;
outputs ... s.` Total timing starts after argument parsing and ends before the
last metadata write and terminal messages. Plotting/GIF time is included under
outputs. See [the architecture notes](docs/runtime-and-parallelisation.md) for
the exact measurement scope and performance interpretation.

Arial is the default requested font throughout static plots and GIFs. If it is
installed, no new option is needed. If it is missing, the terminal announces a
Liberation Sans or DejaVu Sans fallback, and `run.json.plotting` records which
font was actually selected. To require Arial and register a local regular font
file when necessary, use:

```bash
python -m polymer_clustering --config examples/mm32_p3_full.json --strict-font
python -m polymer_clustering --config examples/mm32_p3_full.json --font-path /path/to/Arial.ttf --strict-font
```

The second command uses a placeholder: replace it with the real local font-file
path. The project does not bundle or download fonts. `--strict-font` rejects an
unavailable Arial request before starting a run. `--font-family` permits a
recorded presentation override; font controls never change the simulation.
Saved settings include font options. SVG glyph paths preserve the rendered
appearance; new fonts apply to newly generated figures, not old output files.

Elapsed-time reporting and Arial-first typography from v0.5.1 apply to individual
runs and the new ensemble plots. Parallel execution is now available for whole
independent runs; a single trajectory retains its sequential event algorithm.

## Ten independent replicates per condition (v0.6.0)

The long MM32 comparison is supplied as a batch configuration with ten seeds,
2026–2035, and two conditions: no polymer and 1,000 P3 entities. Validate the
plan, then use the short software diagnostic before launching the long batch:

```bash
python -m polymer_clustering.batch --config examples/mm32_replicates.json --dry-run
python -m polymer_clustering.batch --config examples/replicate_smoke.json --workers 1 --output results/exercise5/batch_check
```

The smoke configuration deliberately accelerates BB170 division and production;
it checks software only. For the actual MM32 experiment, run:

```bash
python -m polymer_clustering.batch --config examples/mm32_replicates.json --workers 1 --output results/exercise5/mm32_replicates
```

Use `--workers 2` if your CPU allocation and memory permit two simultaneous runs.
The default is one. Jobs retain separate logs, timings, arrays and animations;
worker count never selects a seed or alters a trajectory. Ten distinct seeds
are independent within each condition; matching seeds pair the two conditions.
Completed runs can be imported with repeated `--reuse-run` flags, including
compatible v0.5.1 runs, so the existing seed-2026 pair need not be recomputed.
`--resume` retains complete jobs and replays unfinished ones from their seeds.

Only a complete, verified batch receives the ensemble CSV tables and PNG/SVG
comparison. Bands show sample SD across seeds, not confidence intervals. See
the [replicate walkthrough](docs/replicate-walkthrough.md) for the reuse/resume
commands, output interpretation, runtime/resource guidance and scientific limits.

## Begin the affinity controls (v0.7.0)

The five-condition design includes no polymer, P1-like (K_PB=5, K_PS=0), P2-like
(K_PB=0, K_PS=0.1), P3 (5,10), and a P3-without-PB control (0,10). P2's source
signal affinity is 0.1, so it is kept separate from simply disabling P3's PB binding.

Start with the two-seed, 1,000-a.u. pilot:

```bash
python -m polymer_clustering.batch --config examples/mm32_affinity_pilot.json --workers 2 --dry-run
python -m polymer_clustering.batch --config examples/mm32_affinity_pilot.json --workers 2 --output results/exercise5/mm32_affinity_pilot
```

The later `mm32_affinity_full.json` plan contains 50 runs. `--reuse-batch` imports
the twenty matching completed MM32 runs, leaving thirty new simulations. It
copies verified originals, preserves their metadata and never counts a seed
twice. `pairwise_differences.csv` adds descriptive paired contrasts between all
condition pairs when there are three or more conditions. Read the
[affinity walkthrough](docs/affinity-walkthrough.md) for source values, chosen
kinetics, full reuse commands, structural checks and interpretation limits.

## Change controls without editing Python

```bash
python -m polymer_clustering --help
python -m polymer_clustering --seed 2027 --polymers 1000
python -m polymer_clustering --strain BB170 --output ./results/exercise5
python -m polymer_clustering --config examples/mm32_initial.json
```

Command-line values override settings loaded from JSON. JSON values override
defaults. A polymer count is not a physical concentration or a verified conversion
of the paper's `[P]` axis. BB170 starts with zero signals. It produces new signals in full mode; the
earlier diagnostic modes deliberately keep their original fixed populations.

## Run binding

Start with bacterium-signal binding using the published individual rates:

```bash
python -m polymer_clustering --mode binding --output results/exercise5
```

Then exercise all three reaction classes with a documented P3 example:

```bash
python -m polymer_clustering --config examples/mm32_p3_binding.json --output results/exercise5
```

This example has 100 polymer entities and a duration of 100 a.u., chosen for a
quick diagnostic. Its polymer rates are **our explicit choices**, with the
reported P3 ratios `K_PB = 5` and `K_PS = 10`; they are not published individual
kinetic constants. No physical concentration conversion is implied. With polymers
present, all four polymer rates must be supplied, through JSON or CLI flags.

| Pair | Forward flag | Reverse flag | Example forward / reverse |
| --- | --- | --- | --- |
| Bacterium-signal | `--k-bs-on` | `--k-bs-off` | 0.001 / 0.002 (S21) |
| Polymer-bacterium | `--k-pb-on` | `--k-pb-off` | 0.005 / 0.001 (chosen) |
| Polymer-signal | `--k-ps-on` | `--k-ps-off` | 0.01 / 0.001 (chosen) |

For a P1-style binding experiment, override `--k-ps-on 0`. The source P2 ratios
require `--k-pb-on 0 --k-ps-on 0.0001` with the example's reverse rates 0.001.
Disabling only PB retains P3's PS affinity of 10 and is a separate mechanistic
control, rather than the source P2 setting of 0.1. Zero forward rates disable the
corresponding class from an initially unbound state. Positive reverse rates
remain set. A zero reverse rate is allowed for diagnostic experiments
but is an irreversible model change, with no finite association ratio.

`--duration`, `--sample-interval`, `--max-events`, and `--seed` control the
experiment without source edits. Sample times do not change the event sequence.
Increasing duration or polymer count can increase runtime substantially: this is
a readable reference implementation that rebuilds channels and checks invariants
after every chemical event in fixed-position binding mode.

## Make them move

Update the environment from the requirements file, then run either experiment:

```bash
python -m pip install -r requirements-validated.txt
python -m polymer_clustering --mode movement --output results/exercise5
python -m polymer_clustering --config examples/mm32_p3_movement.json --output results/exercise5
```

The first command after installation runs the no-polymer condition. The second
runs the explicit 100-polymer example. Both simulate 100 a.u., with the published
directional movement rates: signals 0.05, polymers 0.01, bacteria 0.0005. Each legal
direction has that rate; do not divide it by four. A cluster of m cells uses
0.0005/m in each legal direction. Movement and binding compete in the same SSA.
The default event cap is now 500000, because signal hops are much more frequent
than binding. Runtime depends on hardware and settings; the included example
took about 12 seconds in the validated environment.

The source leaves transport geometry incomplete. The default obstruction factor
0.1 and the cluster displacement rule are **documented interpretations**, not
measured or uniquely specified source values. Read MOVE-01 through MOVE-07 in
the model notes before interpreting a trajectory scientifically.

| Control | Meaning | Default |
| --- | --- | --- |
| `--d-s`, `--d-p`, `--d-b` | Per-direction signal, polymer, and bacterial rates | 0.05, 0.01, 0.0005 |
| `--obstruction-factor` | Multiplier on small-entity edges touching a bacterial footprint; 0 blocks, 1 removes this obstruction | 0.1 (chosen) |
| `--crowding-lambda` | Centre of the S22 polymer-dependent signal mobility function | 8 |
| `--advection` / `--no-advection` | Enable/disable the documented cluster displacement rule | Enabled |
| `--validate-every` | Full conservation/geometry check interval in events; 1 checks every event | 1000 |
| `--spatial-frames` | Spatial observations for maps and GIF, from 2 to 201 | 21 |
| `--event-log` | `compact` omits signal-hop rows; `all` records every event | compact |

Local non-negativity and receptor capacity are checked on affected sites. Full
state checks run at observations, every `validate_every` events, and at the end.
The cached solver updates affected propensities and uses weighted event trees;
it uses exact SSA for the documented rules, with no fixed integration time step.
Full event logs can be large. Total counts by event type remain in `run.json`
even when signal-hop rows are omitted from the compact CSV. Event numbers in a
compact log therefore have gaps. Change durations or rates through controls;
changing a physical rate is a model experiment and must be reported as such.

## Add delayed expression and predicted light

```bash
python -m polymer_clustering --mode response --output results/exercise5
python -m polymer_clustering --config examples/mm32_p3_response.json --output results/exercise5
```

The first run has no polymers; the second uses the documented 100-polymer example.
Response mode retains binding and movement. Its downstream calculations do not
change the random draws, event history, or spatial arrays. It computes each
cell's response separately; light is then summed/averaged across cells.

| Control | Meaning | Default |
| --- | --- | --- |
| `--k-syn` | Rate at which expression approaches higher sensing | 0.05 |
| `--k-deg` | Positive magnitude of the expression decay rate | 0.01 (sign interpretation) |
| `--response-initial` | Initial expression proxy per cell, from 0 to 60 | 0 (chosen) |
| `--response-midpoint` | Override the light-curve midpoint | 38 MM32; 15 BB170 |
| `--light-scale` | Light at the response midpoint, not the asymptotic maximum | 2500 |
| `--light-map-max` | Fix the light map colourbar across runs | Automatic per-run peak, at least 1 |

At constant sensing S, rising expression E follows `dE/dt = k_syn*(S-E)`.
Falling expression follows `dE/dt = -k_deg*E` until it reaches S, where the
branch changes. These intervals are integrated analytically between occupancy
changes. No Euler time step or per-event expression increment is used. The
supplement prints both a leading minus sign and a negative degradation constant;
RESP-01 records the positive-magnitude interpretation used for decay.

The S22 light formula is retained as printed: `2*Scale / (1 + exp(0.25*(midpoint-E)))`.
Its theoretical asymptote is 5000 for Scale=2500, while at E=midpoint the light is
2500. With E constrained by receptor occupancy to 0-60, that asymptote is not
necessarily attainable. At E=0 the formula still gives nonzero light: about
0.374 per MM32 cell, or 114.887 per BB170 cell. There is no hidden baseline
subtraction. Response mode has no BB170 production; use full mode for that process.
Neither mode alone establishes a completed strain reproduction.

The demonstration's 100-a.u. duration and polymer kinetics remain explicit
choices. A dim response is a valid output. A bright colour in an automatically
scaled map does not imply proximity to the theoretical maximum: read its numeric
colourbar. For comparable visual scales, supply the same `--light-map-max` to
both runs; values above that display limit use the highest colour, while tables
retain their original values. This control changes only presentation.

## Add growth and signal production (v0.5.0)

```bash
python -m polymer_clustering --mode full --output results/exercise5
python -m polymer_clustering --config examples/mm32_p3_full.json --output results/exercise5
python -m polymer_clustering --config examples/bb170_accelerated_growth_demo.json --output results/exercise5
```

The MM32 example retains the 100-a.u. diagnostic duration. With the source's
mean division time of 15000 a.u. and interpreted support 5000-25000 a.u., it will
show no division during that short window. That is expected, not a growth bug.

**The BB170 animation demo accelerates division** to mean 40 a.u., full width 20,
and runs for 150 a.u. It demonstrates births, production, and changing-population
maps quickly. Those growth settings are deliberately different from the paper, and the plots
label that difference.
Production remains 0.005 per cell per a.u., so the short demo may remain dim.
Do not use its timing or light magnitude as a published-figure reproduction.

| Control | Meaning | Default |
| --- | --- | --- |
| `--growth` / `--no-growth` | Enable/disable division in full mode | Enabled |
| `--doubling-time` | Mean of the symmetric triangular cell lifetime | 15000 a.u. |
| `--division-width` | Full support width; zero gives deterministic lifetimes | 20000 a.u. |
| `--population-cap` | Stop division at this population; no deaths or removal | 80 |
| `--division-retry` | Wait before retrying a division with no empty footprint | 100 a.u. (chosen) |
| `--production-rate` | Override signal production per cell; zero disables it | 0 MM32; 0.005 BB170 |
| `--daughter-expression` | Expression proxy assigned to an unbound daughter | 0 (chosen) |

All initial cells are assumed newborn. Each receives a sampled division lifetime;
after a successful birth both the parent and daughter receive fresh lifetimes.
The daughter occupies a nearest available 2 x 2 footprint; tied nearest anchors
are selected uniformly. The parent retains its attachments and expression.
The daughter starts without bound signals/polymers and with its own expression
clock starting at birth. Existing free chemicals in the new footprint remain.
See POP-01 through POP-06 for the exact interpretation choices.

BB170 production is a Poisson event competing with chemistry and transport.
The population rate is cell count times the per-cell rate; each event places one
free signal on one of the producer's four sites uniformly. The checked balance is
`free + bacterium-bound + polymer-bound = initial + produced`. Polymers remain
conserved. No signal degradation, substrate depletion, or cell death is added.

Full mode adds `population_trajectory.csv`, `lineage.csv`, `population.png`, and
`population.svg`. Births, blocked attempts, and production events also appear in
`binding_events.csv`; `bacterium` identifies the parent/producer and `daughter`
identifies a new cell. Initial lineage rows have parent -1 and assumed birth time
zero. Population light sums can rise simply because there are more cells, even
at the light curve's nonzero baseline: compare per-cell means as well as totals.

Frame archives support changing populations without pickle. In
`movement_frames.npz`, unborn position rows are -1 and `bacteria_counts` gives
valid row counts. In `response_frames.npz`, `bacterium` gives the final cell order,
`present` marks existing cells, and unborn response entries are NaN. These NaNs
mean absence, not a zero measurement; CSV response rows list only existing cells.
The maps hide unborn footprints, retain stable cell order, and keep colourbars
fixed across frames.

A source-length MM32 run may require tens of millions of signal-hop events.
The default cap of 500000 can therefore stop well before 30000 a.u.; check
`run.json` for complete status and the actual stop time. Raising the cap is a
runtime choice, not a remedy for unresolved source kinetics or normalisation.

### Placement reproducibility across versions

New CLI runs and new examples use `rng_scheme: independent_v1`. Bacteria, polymer
placement, and signal placement have separate streams. Changing only polymer
count preserves the initial bacteria and free-signal map. This intentionally
changes some seeded initial states from v0.2; it does not reproduce their old
numbers. It also does not guarantee paired event histories after conditions diverge.

Files without `rng_scheme` are treated as v0.2 settings and use `legacy_v02`.
Saved v0.3 settings record the scheme explicitly. Use `--rng-scheme legacy_v02`
for old initial-state replay, or `--rng-scheme independent_v1` to opt an older
configuration into the new placement scheme. The new example files explicitly
select independent streams. Paper-derived settings and algorithm choices remain
separate in the model notes.

Optionally install the package with `python -m pip install -e .`; this also enables
the equivalent `polymer-clustering` command. Module execution above works without
that installation and keeps the first step simple.

## Where outputs go

Each invocation creates a unique child folder beneath `./outputs`, or beneath
the parent specified by `--output`. Paths are relative to the current working
directory. Existing runs are preserved. The program prints the exact folder.

| File | Purpose |
| --- | --- |
| `initial_state.png` | Viewable initial bacterial positions and free-entity maps |
| `initial_state.svg` | Vector version for editing or a report |
| `initial_state.npz` | Occupancy, free signal counts, and free polymer counts |
| `bacteria.csv` | Bacterial identifiers and lower-left lattice coordinates |
| `parameters.json` | Complete effective settings; accepted by `--config` |
| `run.json` | Status, source, interpretation IDs, versions, seed, source hashes, runtime phases, and resolved font |

Binding mode adds:

| File | Purpose |
| --- | --- |
| `binding_trajectory.csv` | Sampled signal pools, attachment counts, bridges, largest cluster |
| `binding_events.csv` | Actual stochastic event times, reaction names, and entity IDs |
| `binding_state.npz` | Final signal pools, polymer identities, and contact arrays |
| `polymer_binding.csv` | Each polymer's anchor, bound signal count, and bacterial contacts |
| `bacteria_binding.csv` | Final bound signals and polymer bonds per bacterium |
| `binding.png`, `binding.svg` | Sampled diagnostic trajectories; no light calculation |

Movement mode retains those common files and adds:

| File | Purpose |
| --- | --- |
| `movement.png`, `movement.svg` | Start/end free-signal maps with bacterial outlines |
| `movement.gif` | Sampled spatial evolution; each frame labels simulation time |
| `movement_frames.npz` | Exact plotted free-signal arrays, times, and bacterial positions |
| `movement_positions.csv` | Bacterial positions at count-observation times |
| `final_bacteria.csv` | Final coordinates, separate from initial `bacteria.csv` |

Response mode adds:

| File | Purpose |
| --- | --- |
| `response_trajectory.csv` | Population sensing/expression and raw light sums/means |
| `cell_response.csv` | Per-cell coordinates, receptor occupancy, expression, and light at sample times |
| `response_state.npz` | Final per-cell response arrays and stop time |
| `response_frames.npz` | Response arrays for each spatial frame |
| `response.png`, `response.svg` | Sensing versus expression, and predicted mean light |
| `response_spatial.png`, `response_spatial.svg` | Final free-signal map beside per-cell predicted light |
| `response.gif` | Paired signal/light maps over time |

These maps do not depict measured photons or experimental microscope images.
They show model outputs. Each coloured 2 x 2 footprint carries one cell's light
value; population light counts each bacterium once. Expression remains associated
with its cell when the cell moves. With no bacteria, sums and means are recorded
as zero as an explicit empty-population convention.

`binding_state.npz` now also stores final bacterial coordinates. In moving event
logs, `x,y` and `to_x,to_y` describe a hop or representative cell translation.
`cluster_members` lists all translated cell IDs; displacement counts describe
free chemicals in the swept region, including chemicals retained at an outer
edge. Chemical event destination fields are blank. Spatial maps show only free
signals, not polymer-bound or receptor-bound populations. Bound populations
remain in the numerical outputs and count trajectories. Animation playback speed
is illustrative and is not a calibration of arbitrary units to hours.

In event tables, `polymer = -1` or `bacterium = -1` means that entity type is not
part of the event. Polymer contact sites are flattened indices: `site = y *
grid_size + x`; contact `-1` means an unoccupied end. A polymer is counted as free
when neither end is bacterium-bound, even if it carries a signal. The bridge count
is the number of polymers attached to two distinct bacteria, not the number of
distinct bacterial pairs. The largest cluster includes isolated cells (size 1).

For example, copy the printed run's `parameters.json` path and supply it to
`--config` to regenerate the numerical initial state in a new output directory.
The same seed, parameters, source, and validated software reproduce the arrays;
timestamps, SVG metadata, and run folder names will differ.

A valid run that later fails leaves `status: failed` and an error in `run.json`
when the folder remains writable. Invalid settings fail before creating a run.
Outputs are ignored by Git by default; select meaningful figures for the report
deliberately rather than committing every exploratory run.

Reaching `--max-events` saves the actual partial history with `status: incomplete`
and exits with code 3. No observations are fabricated beyond the stop time.
Increase the cap and replay the saved settings for a fresh run; this version does
not resume from a saved state. Other exit codes: 0 success, 1 runtime failure,
2 invalid controls. On an unexpected runtime error, metadata is preserved when
possible; partial binding arrays are not checkpointed.

## Understand the code in this order

1. `polymer_clustering/model.py`: configuration, bacterial positions, initialisation.
2. `polymer_clustering/binding.py`: binding state, propensities, events, stochastic loop.
3. `polymer_clustering/output.py` and `binding_output.py`: files, plots, provenance.
4. `polymer_clustering/cli.py`: user controls and the run lifecycle.
5. `polymer_clustering/movement.py`, `rates.py`, and `movement_output.py`: coupled
   transport rules, cached rate selection, spatial views, and animation.
6. `polymer_clustering/response.py` and `response_output.py`: continuous response,
   population/per-cell tables, light plots, and paired maps.
7. `polymer_clustering/population.py` and `population_output.py`: division
   appointments, daughter placement, production, lineage, and population balance.
8. `polymer_clustering/timing.py` and `plotting.py`: real execution measurements,
   scoped fonts, and typography provenance across all graphical outputs.
9. `polymer_clustering/batch.py` and `batch_output.py`: explicitly seeded process
   scheduling, retained attempts, provenance checks, and across-run statistics.
10. `polymer_clustering/progress.py`: read-only time/event reporting and terminal
    or worker-status output, kept outside model parameters.

`occupancy[y, x]` is -1 for an empty bacterial site, otherwise a bacterial ID.
Every bacterium occupies exactly four sites. Independent count arrays allow
signals and polymers to share sites with bacteria. These arrays are not maps of
measured signal concentration; they are counts in the reconstructed initial state.

## Validate

```bash
python -m unittest discover -s tests -v
python -m compileall -q polymer_clustering tests
```

Checks cover overlap, boundaries, counts, reproducibility, invalid settings,
JSON/CLI precedence, failed runs, and preservation of existing output folders.
Binding checks add conservation, capacities, bridge identity, crowding, event
timing, reversible analytical behaviour, sampled observations, and partial runs.
Movement checks add random-stream independence, boundary/overlap restrictions,
carrier transport, collective motion, displacement, cache equivalence to full
recalculation, the free-diffusion law, animation output, and legacy replay.
Response checks add analytical rise/decay, branch crossing, light midpoint and
asymptotes, nonzero baseline, per-cell aggregation, passive coupling, and
observation-independent results. Population checks add lifetime statistics,
scheduled births, parent attachment retention, daughter response clocks, production
statistics/conservation, population caps, variable-size frames, and replay. See [validation](docs/validation.md) for scope.

## Model audit and next steps

Read [the model notes](docs/model-notes.md) for source pages, explicit assumptions,
and unresolved equations. Record changes in [CHANGELOG.md](CHANGELOG.md) and
AI use in [the assistance record](docs/ai-assistance.md).

Read [the binding walkthrough](docs/binding-walkthrough.md) to follow one event.
Read [the movement walkthrough](docs/movement-walkthrough.md) and
[the response walkthrough](docs/response-walkthrough.md) for the new components.
Read [the population walkthrough](docs/population-walkthrough.md) for the new
layer. The ten-seed MM32 affinity comparison, two-seed crowding pilot and focused
concentration pilot have been reviewed. The concentration-dependent switch from
early suppression to later enhancement motivates further analysis of existing
outputs. The current release improves the solver before further long collection.
BB170 initial-expression diagnostics, selected ten-seed extensions and
source-figure checks remain.

No publication PDFs are redistributed. No software licence has been selected:
check the existing assignment repository's licence before public redistribution.

