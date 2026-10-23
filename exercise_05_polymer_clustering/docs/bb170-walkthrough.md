# BB170 data collection for Exercise 5 Part 3

Version 0.11.1 prepares BB170 collection on the student's Linux laptop. The
optimised solver, progress, runtime records, independent-run parallelism,
Arial-first figures and stable-scale animations are retained. All assignment
data for this stage should be generated locally with the recorded configuration
and seeds. The first pair is a reusable beginning, not the final average.

## Source conditions and complete collection

The supplied supplement specifies BB170 initial signals zero (S21), light
midpoint 15 (S22), and signal production 0.005 per cell per a.u. in our documented
Poisson interpretation (S23). S26 specifies polymer values 0-7000 in steps of
1000, and time 0-30000 sampled every 100. S40-S42 use P3 affinities K_PB=5 and
K_PS=10 and change initial density or mean doubling time.

| Plan | Initial bacteria | Mean doubling time | Polymer values | Seeds | Runs |
| --- | ---: | ---: | --- | --- | ---: |
| Standard first pair | 25 | 15000 | 0 and 1000 | 2026 | 2 |
| Standard full, S40 | 25 | 15000 | 0,1000,...,7000 | 2026-2035 | 80 |
| Lower density full, S41 | 2 | 15000 | 0,1000,...,7000 | 2026-2035 | 80 |
| Faster growth full, S42 | 25 | 10000 | 0,1000,...,7000 | 2026-2035 | 80 |

The three full plans contain 240 distinct requests. The first pair is an exact
subset of the standard plan and can be imported after completion; it does not
add two more independent replicates. Separate batches keep different density
and growth regimes separate in tables, plots and averaging.

Every plan uses the full 30000-a.u. duration, 100-a.u. sampling and 61 spatial
frames. The grid remains 50 x 50, the population cap 80, and the PS cutoff ten.
PB on/off rates remain 0.005/0.001 and PS on/off 0.01/0.001, giving the stated
ratios. These individual polymer rates are our documented choices, since the
source does not provide both rates. Signal production is explicitly 0.005 and
the response midpoint explicitly 15 in the JSON, rather than inferred silently
from the strain name.

S26's text conflicts with some 70000-polymer legends and 3000-time captions.
These configurations follow its explicit 0-7000 and 0-30000 statements. Counts
are model entities without an established conversion to experimental mg/ml.

## Lifetime boundary required by the faster-growth condition

Retaining the S23 full width 20000 when S42 changes the mean to 10000 gives a
symmetric triangular distribution on [0,20000]. S42 does not restate the width;
keeping it fixed is an explicit interpretation of the unchanged S23 setting,
not confirmation of the authors' implementation.

The earlier validator rejected equality between the mean and half-width.
Version 0.11.1 accepts this nonnegative support. A continuous triangular
distribution has no probability mass at its zero endpoint. A finite-precision
zero draw is resampled, ensuring every scheduled lifetime is positive.
Negative support and a zero mean remain invalid.

Previously valid positive-support configurations use exactly the same NumPy
draws and random-stream states as before. The default MM32 and standard BB170
support remains [5000,25000]. Existing optimisations and other numerical modules
are unchanged. Reviewed source fingerprints preserve compatible completed-run
reuse only for the old positive-support parameter domain.

## Initial expression remains an explicit assumption

Initial-cell and daughter expression remain zero. S21-S22 do not specify their
initial values, so this choice remains recorded rather than selected to fit a
curve. The source light transform, including its nonzero baseline, is retained.
This uncertainty belongs in the later Part 4 comparison; these plans do not
claim exact agreement with the paper.

## 1. Start the first full-duration pair on the laptop

Extract v0.11.1 to a new folder, enter its polymer-clustering directory and
activate the existing validated environment. Dependencies have not changed.
Check the source version and preview:

```bash
python -m polymer_clustering --version
python -m polymer_clustering.batch --config examples/bb170_standard_pair.json --workers 2 --dry-run
```

Expect version 0.11.1 and **2 pending runs**: no polymer and P3 at 1000 entities,
both seed 2026. The preview starts no solver and creates no run folders.

Then run:

```bash
python -m polymer_clustering.batch --config examples/bb170_standard_pair.json --workers 2 --progress on --output results/exercise5/bb170_standard_pair
```

Keep the unique batch folder printed at launch. Each trajectory runs on the
laptop and records the source fingerprints, NumPy, seed, initial conditions,
production, growth, output font and timings. The batch is complete only when
both trajectories and their outputs have passed validation.

## 2. Review the first pair

Send these five files from the printed batch folder:

- batch.json
- ensemble.png
- ensemble_trajectory.csv
- replicate_summary.csv
- paired_differences.csv

An intact small archive is convenient. From the project directory, replace the
placeholder with the actual printed batch folder:

```bash
PAIR_BATCH="/actual/printed/completed/pair/batch"
tar -czf bb170_standard_pair_summary.tar.gz -C "$PAIR_BATCH" batch.json ensemble.png ensemble_trajectory.csv replicate_summary.csv paired_differences.csv
```

Keep the individual trajectory CSVs, NPZ states and GIFs locally. We can select
them if a check needs them. This first pair has one seed per condition: its mean
is that trajectory, SD/SEM are undefined (blank in the tables), and no uncertainty
band is drawn. Ten independent seeds are still required for the final average.

Review completion and conservation first, then signal production, receptor
sensing, expression, per-cell light, population and clustering. BB170 has
`free + bacteria-bound + polymer-bound signals = produced signals` at every
observation, because it starts with zero and has no specified signal loss.
Population cannot exceed 80; receptor and polymer capacities are also checked
by the existing validator. Part 4 will compare the collected curves with the
publication after Part 3 collection is complete.

## 3. Extend to the standard ten-seed loading sweep

Use the actual complete first-pair folder. Preview the strict reuse:

```bash
PAIR_BATCH="/actual/printed/completed/pair/batch"
python -m polymer_clustering.batch --config examples/bb170_standard_full.json --reuse-batch "$PAIR_BATCH" --workers 6 --dry-run
```

Expect **80 jobs: 2 reused and 78 pending**. Parameters, complete outputs,
numerical source identities, NumPy and observation grids must match. The
controller copies verified originals and never counts a seed twice. No MM32
or workspace-generated BB170 run is imported into this laptop collection.

Then launch:

```bash
python -m polymer_clustering.batch --config examples/bb170_standard_full.json --reuse-batch "$PAIR_BATCH" --workers 6 --progress on --output results/exercise5/bb170_standard_full
```

Six workers follow the student's previously used setting for a computer now
dedicated to simulation. Worker count affects scheduling, not scientific
parameters or seeds. It is not a promise of linear speedup or a completion
deadline. Existing per-run timings and progress show actual performance.

Only all ten complete, verified seeds per loading receive the final ensemble
average. Pointwise bands show sample SD across seeds, not a confidence interval.
The complete loading sweep adds pairwise_differences.csv to the review files.

## 4. Collect the other reported BB170 regimes

After reviewing the standard study, preview and run lower density:

```bash
python -m polymer_clustering.batch --config examples/bb170_low_density_full.json --workers 6 --dry-run
python -m polymer_clustering.batch --config examples/bb170_low_density_full.json --workers 6 --progress on --output results/exercise5/bb170_low_density_full
```

Then faster growth:

```bash
python -m polymer_clustering.batch --config examples/bb170_faster_growth_full.json --workers 6 --dry-run
python -m polymer_clustering.batch --config examples/bb170_faster_growth_full.json --workers 6 --progress on --output results/exercise5/bb170_faster_growth_full
```

Each full regime has 80 pending runs when started fresh. Use one batch at a time
so the total worker count remains six. A full sweep can be reviewed using the
five files listed above plus pairwise_differences.csv; large individual outputs
need not all be uploaded.

## Interruptions and completion

Resume an unfinished batch with the same unchanged package and environment,
using its printed batch folder:

```bash
python -m polymer_clustering.batch --resume "/actual/printed/batch/folder" --workers 6 --progress on
```

Completed jobs are retained; unfinished trajectories replay from their seeds.
A simulation reaching its event budget is incomplete and is not averaged.
Never edit metadata to turn an incomplete trajectory into a complete one.

Part 3 BB170 collection is complete for these targets when all three full
regimes have eight loadings and ten verified independent seeds per loading.
Software diagnostics, the first illustrative pair and any intermediate subset
do not replace that averaging. Keep the exact package, configurations, output
folders and logs for subsequent Part 4 comparisons.

