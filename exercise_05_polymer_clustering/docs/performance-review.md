# Solver performance review - 22 October 2026

Version 0.11.0 addresses an algorithmic bottleneck found after the student's
7000-polymer runs took about fourteen hours each. The improvement applies to
cache bookkeeping. Equations, rates, geometry, source assumptions, seeds,
SSA event ordering, sampling, validation frequency and animations are retained.

## Evidence from the completed concentration batch

The uploaded v0.10.0 records contain four reused and four new runs. The two
3000-polymer runs consumed 6.59-6.67 CPU hours for about 7.6 million events
each. The two 7000-polymer runs consumed 14.03-14.20 CPU hours for about
6.4 million events each. Their average CPU cost was approximately 3.15 versus
7.99 milliseconds/event. More events cannot explain this slowdown.

These are mixed-workload averages, not per-event-kind measurements. They do not
by themselves identify the slow function. Individual CPU time is close to the
recorded elapsed time, which supports investigating CPU work rather than
assuming an I/O or worker-count problem. Four pending trajectories allowed at
most four active workers initially; after two finished, only two remained.

## Confirmed cause

`MovementEngine.event()` called `retether()` over every polymer and then
`rebuild()` after every PB bind/unbind. Rebuild regenerated neighbour/contact
information, local chemical channels, diffusion rates, all weighted trees and
cluster movement rates across the entire state. This happened even when one
contact changed and all bacterial positions were unchanged.

A deterministic cProfile check of the first 5000 events on the 50 x 50,
25-bacterium MM32 setup found 126 global rebuilds at 7000 polymers, including
initialisation. About 96.5% of the profiled solver time was inside rebuild.
The profile counted 905643 calls to polymer_reactions for those 5000 events.
The optimised check needs one initial geometry rebuild and 31145 calls to that
function for the same event sequence. Imports and plotting are excluded.

## Changes

1. Index every structurally possible PB target, including targets at capacity
   whose positive channel is currently absent. A changed target refreshes its
   dependent carriers, including singly bound carriers on neighbouring hosts.
2. Retether only the changed carrier. If its anchor moves, update both local
   bins, PS accessibility and free-signal mobility. Polymer hops maintain the
   same dependency index and mark the insertion-order cache dirty.
3. Recompute connectivity only when a second end changes. Recalculate cluster
   motion rates when the connected components change. Actual translations and
   divisions still rebuild the full geometry.
4. Preserve v0.10.0 floating-point sums after PB events. Cached leaf rates are
   inserted into zero Fenwick trees in the old site/polymer order, using
   unbuffered NumPy addition and bounded cached integer paths. A differently
   ordered bulk sum could otherwise alter seeded event times or selection.
5. Skip unchanged-rate ancestor walks. Mask bridge rows in one NumPy operation
   before the existing ordered connectivity traversal.

The dependency cache is per process. Integer Fenwick paths use a bounded
eight-size cache, rather than retaining an unbounded history of simulation
states. There is no shared mutable state between replicates and no new dependency.

## Measurements

Same Linux host, Python 3.12.14, NumPy 2.3.5 and the pinned dependencies.
Each table uses the median current-process CPU time across three **unprofiled**
repetitions. Profiling is a separate repeat, so instrumentation overhead does
not enter these timings. Times omit imports, output files and rendering.
The baseline and candidate were tested on the same host; ordinary frequency/load
variation remains possible. These measurements are reproducible diagnostics,
not hardware-independent benchmarks.

First 5000 exact SSA events, source-default initial grid/bacteria/signals/P3
kinetics, seed 2026, full response/population mode, all events logged:

| Polymers | v0.10.0 CPU seconds | v0.11.0 CPU seconds | Speedup |
| ---: | ---: | ---: | ---: |
| 1000 | 1.949 | 0.283 | 6.9x |
| 3000 | 7.328 | 0.519 | 14.1x |
| 7000 | 21.522 | 0.849 | 25.3x |

Each trajectory deliberately stops at its 5000-event cap. It is incomplete as a
scientific run and samples the early state, rather than the full 30000 a.u.
For every count, the old/new result digest is identical over **all result fields**,
including every logged event/time, state-array byte, observation, per-cell
response, frame, lineage entry and event count.

A second test uses an explicitly **constructed valid crowded state**, with
79 bacteria in four separated packs of at most twenty cells, many bridges,
occupied PS slots and 10000 conserved signals. It is not an uploaded endpoint
or a predicted scientific result. Five valid bind/release cycles (ten PB events)
exercise the dense cache path and restore the physical state:

| Polymers | v0.10.0 CPU seconds | v0.11.0 CPU seconds | Speedup |
| ---: | ---: | ---: | ---: |
| 1000 | 0.567 | 0.0124 | 45.8x |
| 3000 | 1.047 | 0.0307 | 34.1x |
| 7000 | 1.995 | 0.0753 | 26.5x |

Those old/new digests include identical final arrays, component identities,
rate leaves and the internal weighted partial sums. Raw timings, software
fingerprints, digests and profile summaries are in `benchmarks/results/`.

The full fourteen-hour trajectories have **not** been rerun here. Late density,
bridging, moves, division, logging and rendering change the work mix. Do not
divide the old runtime by a short-test speedup to predict completion time.
The supplied new summaries do not include their full endpoint arrays, so the
crowded fixture is explicitly labelled rather than substituted for real data.

## Correctness and reuse

All 140 tests pass. New checks compare dependent caches and exact weighted trees
with independent full reconstruction after capacity release/refill, bridge
formation/removal, reanchoring across both sides of the PS cutoff, redundant
bridges, polymer hops and insertion-order changes. Randomized rates at several
sizes/magnitudes check exact old insertion arithmetic. Mixed growth, production,
transport and response runs also replay a full-rebuild oracle exactly.

Complete v0.10.0 CLI references in binding, full MM32, cutoff-off and accelerated
BB170 modes were replayed. Numerical CSVs are byte-identical and every saved
array byte matches (11 arrays in binding mode, 27 per coupled mode). These small
checks include division/production; they are not assignment-scale predictions.

Among the six numerical modules, only movement.py and rates.py change. Model,
binding, response, population, pinned requirements and dependency ranges remain
byte-identical to v0.10.0. Two exact fingerprint bridges permit verified
completed-run import. The earlier binding/default-ten bridges remain restricted
to their reviewed hashes and conditions. Unknown edits and NumPy differences
are rejected; completion, parameters, grids and balances are still checked.

The uploaded eight-job concentration manifest passes the source compatibility
check. Its complete arrays are on the student's computer and remain subject to
the existing local import validation. There is no need to repeat the completed
concentration pilot. Resume of unfinished work requires its original complete
package/environment; cross-version reuse of finished runs is separate.

## Reproduce the review on Linux

Extract v0.11.0 to a new folder and activate the existing environment. Keep the
old extracted v0.10.0 folder for comparison. No dependency update is needed.
From the new project folder, replace the old-folder placeholder:

```bash
OLD_PACKAGE="/absolute/path/to/v0.10.0/polymer-clustering"

python benchmarks/profile_solver.py \
  --source-root "$OLD_PACKAGE" \
  --case trajectory --polymers 1000 3000 7000 \
  --events 5000 --repeats 3 --profile \
  --output results/exercise5/performance/v010

python benchmarks/profile_solver.py \
  --case trajectory --polymers 1000 3000 7000 \
  --events 5000 --repeats 3 --profile \
  --output results/exercise5/performance/v011
```

The output folders are unique. Compare `result_sha256` for each polymer count;
it must match across versions. Compare median CPU seconds with the same inputs
and environment. Use `--case crowded-pb --cycles 5 --repeats 3 --profile` in both
commands for the dense-state test. `pb-cycle` supplies the simpler fresh-state
cache diagnostic. Every profile is saved separately as a `.prof` file for
`pstats`; the companion JSON includes readable top-function summaries.

## Remaining limits and next targets

Exact SSA trajectories remain sequential: the next event depends on the state
and rates after the previous one. Worker parallelism applies to independent
trajectories, not the internal events of one trajectory. More workers cannot
accelerate two remaining jobs by themselves and can compete for CPU/memory.
The current implementation respects Linux CPU affinity/SLURM allocation and
limits nested numerical-library threads inside each worker.

Dense sites still require refreshing multiple local polymer channels after a
free-signal change. Bridge-component checks, genuine geometry changes, weighted
selection, state validation and per-cell observations also cost CPU. Profile
representative workloads again before extending the next optimisation. Rendering
and large event logs remain separate costs. They were not removed to improve
the solver benchmark, and no faster approximate stochastic method was introduced.

Future larger sweeps should begin with a bounded diagnostic and a representative
pilot. Keep numerical signatures, solver/output phase times and per-run CPU
measurements. The current result establishes a substantial avoidable cost,
without promising a particular complete-run time on another computer.

Software references: Python's [cProfile guidance](https://docs.python.org/3.12/library/profile.html)
and NumPy's [unbuffered ufunc.at](https://numpy.org/doc/2.3/reference/generated/numpy.ufunc.at.html).
These explain tooling/arithmetic, not biological rules.
