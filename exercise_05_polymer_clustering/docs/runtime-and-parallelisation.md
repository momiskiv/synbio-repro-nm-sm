# Execution time, typography, and parallelisation

Version 0.5.1 adds execution measurements and a shared plotting policy. It does
not change reaction propensities, event selection, division clocks, or expression.
Version 0.6.0 adds the batch controller described below. The numerical solver
inside each run remains sequential. These are software architecture changes,
separate from the source-model interpretations in model-notes.md.

## Timing contract

`timing.py` uses Python's monotonic `time.perf_counter` for elapsed wall time.
This measures how long the user waits, including time spent waiting for the CPU
or files. `time.process_time` separately measures CPU work by the current process.
Neither is a conversion of the simulation's arbitrary time units into seconds.

The timer starts after command-line parsing. Its final snapshot is taken after
all requested arrays, tables, figures, and GIFs have been written, immediately
before the final `run.json` write. Interpreter imports, argument parsing, the
last metadata write, and terminal printing are excluded. UTC timestamps remain
useful for chronology; clock-adjustable UTC differences are not used as timers.

Every started run records these `run.json` fields on completion, an event-limit
stop, or a caught runtime failure, provided the metadata file can be written:

| Field under `timing` | Meaning |
| --- | --- |
| `total_wall_seconds` | Elapsed time over the scope above |
| `initialisation_wall_seconds` | Generate the numerical initial state |
| `simulation_wall_seconds` | Binding solver or coupled solver, including its response calculations and observations |
| `output_wall_seconds` | Initial and final file output, plotting, and GIF encoding; includes initial metadata writes |
| `other_wall_seconds` | Remaining setup, configuration, font resolution, and metadata preparation |
| `cpu_seconds` | CPU time consumed by this process |
| `simulation_events_per_second` | Actual executed events divided by simulation wall time, or null without a completed solver result |

Phases do not overlap. Repeated output phases accumulate. A phase's `finally`
block captures partial elapsed work if it fails. An incomplete run's throughput
uses its actual executed event count, not the requested event cap. CPU time does
not include separate processes. The v0.6.0 batch controller preserves each
run's measurements and reports its own wall time separately. Invalid controls
do not start a run or create outputs.

Terminal output gives total wall time, simulation time, and output time. Longer
totals also show hours:minutes:seconds. Initialisation and the finer breakdown
remain in metadata. Timing is diagnostic, not a hardware-independent benchmark;
CPU load, storage, font caches, library versions, and plotting choices matter.

## One typography policy

`plotting.py` selects Arial by default for titles, labels, legends, ticks,
colourbars, and GIF text. All six plotting entry points use its scoped decorator;
the CLI supplies one resolved style for the whole run. Matplotlib settings are
restored afterwards, including if rendering fails, so importing this project
into another program does not permanently change its plot defaults.

The project uses an installed font or registers the local file supplied through
`--font-path`. The file's internal family must match `--font-family`. Fonts are
not downloaded or redistributed. If Arial is unavailable, the program announces
and records Liberation Sans, or bundled DejaVu Sans, as a fallback. It never
labels that fallback as Arial. `--strict-font` rejects a missing requested font
before creating output files. An explicit font file is useful on machines that
cannot install fonts system-wide.

Settings `font_family`, `font_path`, and `strict_font` are saved in
`parameters.json`. `run.json.plotting` records requested and resolved family,
actual font path/hash, whether fallback occurred, and any registered file.
Relative font paths resolve from the working directory. Replaying a file-based
font request on another machine requires making that font file available there.
SVGs retain glyph paths for a consistent appearance when viewed elsewhere;
PNG/GIF lettering is rasterised using the resolved font. Earlier saved figures
are not retroactively restyled; numerical outputs remain suitable for comparison.

## Parallelisation decision

**Implemented in v0.6.0: parallelise independent runs in a batch runner.**
Each individual exact-SSA trajectory remains sequential. Each event
changes counts and rates used by the next event. Simply distributing those
updates across threads would not preserve the existing event algorithm.

The assignment's ten independent replicates and polymer-condition sweeps are a
useful unit of concurrent work. The controller launches the existing guarded
single-run CLI in fresh subprocesses using the current Python interpreter. This
avoids the GIL restriction on Python-heavy work, shared Matplotlib/RNG state,
and transferring large frame/event arrays through a process pool. Its parent
only schedules, polls, validates and summarises small trajectory tables. It does
not run solver events. `--workers` is explicit and defaults to one.

The supplied seed-2026 MM32 pilot illustrates why profiling comes first:

| Condition, t=1000 a.u. | All executed events | Signal hops | Previous recorded elapsed UTC interval |
| --- | ---: | ---: | ---: |
| No polymer | 1,751,510 | 1,747,888 | 84.219 s |
| 1000 P3 entities | 1,003,073 | 950,701 | 76.734 s |

These are historical v0.5.0 start-to-finish intervals from the student's run
records, not the new phase measurements or a controlled speed benchmark.
Polymer binding reduced free signals and therefore hop events, while adding
polymer reactions and movement work. Runtime is not necessarily monotonic in
polymer count. Longer durations, more mobile signals, geometry, attachment
changes, recording frequency, and large logs/animations can all affect cost.

The batch controller follows these requirements:

1. Give each job an explicit seed and full settings; independent jobs do not share
   random generators. Job ordering or worker count must not alter a seed's result.
2. Keep each run's unique output directory, status, timing, and font provenance.
   Return compact summaries rather than transferring large frame/event arrays
   between processes. Plot within each isolated worker or render separately.
3. Limit workers to allocated CPU and memory resources. On an HPC system, use the
   scheduler allocation, not the machine's total advertised core count. Additional
   workers duplicate state, frame buffers, and logs and may compete for disk I/O.
4. Keep a serial execution option and test the same seed list with one and several
   workers. Numerical arrays/event counts must agree per job; timings will differ.
5. Benchmark equal workloads with identical output settings. Report total batch
   wall time, per-run simulation/output timings and incomplete/failed jobs. Memory
   is not measured automatically; assess it separately before increasing workers.
   Parallel runs do not guarantee linear speed-up.
6. Exercise actual separate interpreters. Children run the single-run CLI only;
   they never launch batches recursively. Linux subprocess execution is tested.
   The controller requires POSIX file locking; Windows is not supported, and
   macOS execution has not been validated in this environment.

See [the replicate walkthrough](replicate-walkthrough.md) for the interface,
source/NumPy checks on imports, retained attempts, resource limits and timing
scope. A batch can replay unfinished jobs; it cannot restart an SSA from an
internal checkpoint. Existing single-run commands retain their behaviour.

The student's completed affinity workload used six workers on twelve available
logical CPUs but made concurrent writing uncomfortable. Short CPU/RAM snapshots
did not establish sustained desktop responsiveness. For this machine, three or
four workers are a more conservative starting point during interactive work;
six can be reserved for unattended batches. Compare resource use over a full
workload, and do not infer runtime speed-up from different condition sets.

The v0.8.0 spatial diagnostic is post-processing, not a worker batch. It reads
saved final arrays, verifies their provenance/counts, and renders a small number
of figures. Its terminal reports elapsed time; spatial_report.json records wall
and CPU seconds with scope. No simulation is executed.

## Progress reporting - version 0.10.0

Both execution CLIs accept `--progress auto|on|off`. Auto redraws one line only
in a terminal. On also writes plain periodic lines when redirected, without
ANSI or carriage-return sequences. Off hides progress; normal lifecycle/timing
messages remain. No new dependency is needed. Display mode is execution metadata,
never a model parameter or a reason to reject a matching completed trajectory.

Single runs show initialisation, initial outputs, simulation, final outputs and
status. Simulation percentage is actual model time divided by requested duration.
Optional read-only callbacks observe time/event count at start, every 1000 events
and stop. They add no simulation observations, rate updates, random draws or
division appointments. Display/file writes are throttled using monotonic real
time. With display off and no parent status file, the solver has no callback.
Counts include all signal hops, regardless of compact/all event logging.

Displayed batches use small atomic per-attempt `progress.json` worker files,
updated at most once per real second plus phase changes. Workers print no live
bars into console logs. The controller shows latest phases and verified complete
jobs over all jobs, including imports. There is no remaining-time prediction.
Progress writes are best-effort and cannot change the numerical trajectory.
Status files are optional, not restart checkpoints or final scientific metadata.

Simulation 100% precedes rendering. Jobs count complete only after successful
exit and trajectory validation; ensemble rendering is announced separately.
Event-limited, failed and interrupted jobs remain labelled and incomplete.
Ctrl-C terminates batch children, retaining complete jobs and prior attempts.
Standalone Ctrl-C also records interruption status/timing when possible.

Only two of six numerical files change, for progress hooks/comments. Exact
fingerprint bridges accept reviewed binding/movement predecessors; arbitrary
edits are rejected. Legacy fixed-cutoff binding still requires ten; the v0.9.0
predecessor supports its saved threshold, including off. NumPy and the other
four numerical hashes must match. Copied original metadata stays intact.
Resume retains its original full-package/environment requirement; completed-run
cross-version reuse is separate.

The concentration pilot recommends two workers initially while writing. This
is a conservative resource choice, not a measured speed-up. Higher loading has
no established runtime multiplier. Retain timings and assess responsiveness.

## Profiled solver optimisation - version 0.11.0

The concentration workload identified a substantial algorithmic cost: PB bond
changes rebuilt every geometry/rate cache despite fixed bacterial positions.
The new structural PB-target index refreshes dependent polymers, local PS/mobility
rates and changed connected components. It includes zero-rate capacity targets
and maintains dependencies after carrier movement. Genuine translations and
division retain their full geometry rebuild.

Weighted sums are recreated using the old insertion order in unbuffered NumPy
operations, rather than permitting a new floating-point accumulation history.
This preserves the previous exact seeded event/time sequence in the verification
cases. Integer update paths have a bounded per-process cache. The existing
validation frequency, logging modes, plotting, progress and worker isolation stay
available. Four numerical modules and all dependencies remain byte-identical.

The [performance review](performance-review.md) records the measured bottleneck,
short and crowded benchmarks, exact result digests, Linux reproduction commands,
remaining costs and limitations. Exact SSA events remain sequential; independent
replicates are the supported unit of process parallelism. Complete old data can
be imported through exact reviewed fingerprints. Pending batches must resume
with their original package/environment.

For Linux timing interpretation, the current monotonic perf-counter clock
excludes system suspend. Ctrl-Z stopping a process while the OS remains awake
continues to accumulate elapsed time; process CPU time does not. UTC differences
also include sleep but may reflect clock adjustments. Preserve both measurements
instead of labelling all elapsed time as CPU processing.

## Implementation references

- Python 3.12 timing: https://docs.python.org/3.12/library/time.html
- Python 3.12 process pools: https://docs.python.org/3.12/library/concurrent.futures.html
- Matplotlib 3.10.8 font discovery/registration: https://matplotlib.org/3.10.8/api/font_manager_api.html

These document software behaviour; they do not supply biological model rules.
