# Changelog

## 0.11.1 - 22 October 2026 - BB170 collection and lifetime boundary

- Added a reusable full-duration BB170 first pair and three full 80-run plans:
  standard, lower initial density and faster growth, each with P3 values
  0-7000 by 1000 and ten seeds (2026-2035).
- Allowed a continuous triangular lifetime support starting at zero, needed
  for mean 10000 with retained width 20000. Finite-precision zero samples are
  redrawn; a zero mean or negative support remains invalid.
- Kept old positive-support sample values and RNG states exactly unchanged.
  Added a narrowly restricted, exact population-source reuse bridge.
- Added boundary/distribution, coupling, reproducibility, compatibility and
  collection-plan checks. Dependencies and existing solver optimisations remain.
- Labelled the BB170 faster-growth condition accurately in individual and batch
  figures, and made singleton ensemble titles explicit about missing SD.
  All eight loading conditions now have distinct colours and explicit labels.
- Documented laptop collection, full-duration reuse, source ambiguities and
  separate Part 3 collection versus later Part 4 comparison.

## 0.11.0 - 22 October 2026 - Profiled solver optimisation

- Replaced whole-state PB-event rebuilds with structural target dependencies,
  affected-site chemistry/mobility updates and conditional component updates.
  Capacity-zero channels, polymer hops and reanchoring remain accounted for.
- Preserved the old Fenwick insertion arithmetic through ordered unbuffered
  bulk updates; skipped unchanged-rate walks and vectorised bridge filtering.
- Added bounded fresh/crowded profiling diagnostics, exact result digests, raw
  benchmark records and a performance review. Repeated short trajectory tests
  measured 6.9x/14.1x/25.3x speedups at 1000/3000/7000 polymers on this Linux host.
- Added eight optimisation/compatibility tests; retained biological rules,
  observations, validation intervals, output quality and dependencies.
- Added reviewed movement/rates fingerprint bridges for complete-run reuse.
  Resume still requires the original full package/environment. No complete
  fourteen-hour trajectory was rerun and no full-run speedup is promised.

## 0.10.0 - 21 October 2026 - Progress and concentration pilot

- Added auto/on/off terminal progress, plain redirected output, simulation/output
  phases, verified batch-job counts and optional atomic per-worker status.
- Added read-only solver hooks without changes to event scheduling, random draws,
  biological rules or dependencies; retained failure/interruption/event limits.
- Added an eight-trajectory MM32 pilot at 0/1000/3000/7000 entities, seeds
  2026/2027 to 30000 a.u.; four compatible prior runs leave four new simulations.
  Added concentration labels, source notes and a collection walkthrough.
- Added exact predecessor fingerprint bridges for completed-run reuse; legacy
  fixed-cutoff runs remain restricted to ten. Resume stays package-specific.
- Recorded standalone interruption status/timing and added seeded numerical
  equivalence, progress and configuration checks. No long scientific pilot was
  run for software verification.

## 0.9.0 - 21 October 2026 - Matched PS crowding sensitivity

- Parameterised PS-site occlusion with a positive threshold or explicit off/null;
  retained default ten, PB kinetics, PS release and transport rules.
- Added a four-run, two-seed MM32 pilot to 10000 a.u., changing only PS occlusion
  between conditions. Outputs label altered-rule sensitivity.
- Made spatial diagnostics respect each run's cutoff and retain fixed P >= 10
  density metrics, with active-cutoff labels on shared-scale maps.
- Added a narrowly fingerprinted default-ten bridge for completed legacy runs.
  Five numerical modules and dependencies remain unchanged; default CSV/array
  replay was verified against v0.8.0 runs.

## 0.8.0 - 21 October 2026 - Endpoint spatial diagnostics

- Added an independent read-only spatial-analysis CLI for complete batches,
  using saved polymer identities, positions, contacts and signal counts. No
  simulation replay or new multi-hour ensemble is needed.
- Added endpoint density, PS crowding/accessibility, localisation, capacity
  occupancy and instantaneous-rate diagnostics; same-state rates without the
  cutoff are explicitly hypothetical, not intervention trajectories.
- Added validated per-run/site data, equal-weight ensemble and paired summaries,
  Arial-first shared-scale maps, elapsed time, input hashes and a small review ZIP.
- Documented that existing frame arrays cannot supply polymer time courses.
  The six numerical solver modules, source rates and dependencies are unchanged.

## 0.7.0 - 20 October 2026 - Affinity controls

- Added matched MM32 pilot/full configurations for P1-like (5,0), P2-like (0,0.1),
  P3 (5,10), P3-without-PB (0,10), and no-polymer conditions. Documented source
  affinities separately from the chosen individual kinetics and experiment scope.
- Added `--reuse-batch` to import all matching complete trajectories from a prior
  complete batch, with existing per-run checks, duplicate detection, preserved
  original metadata, path containment, and dry-run reuse/pending counts.
- Added paired contrasts between every condition pair for multi-condition batches,
  readable affinity labels and consistent colours. Existing baseline contrasts,
  SD/SEM definitions, timing and Arial policy remain available.
- Added structural control and reuse/contrast checks, a staged walkthrough and
  an ignore rule for the default results directory. No model rates/rules were
  changed in the six numerical source modules; dependencies remain unchanged.

## 0.6.0 - 20 October 2026 - Independent replicate batches

- Added explicit condition/seed batch plans and optional bounded child processes;
  every numerical SSA trajectory retains its serial algorithm and random streams.
- Added dry-run validation, per-attempt console logs, batch progress/timing,
  Ctrl-C child cleanup, retained attempts and replay of unfinished batch jobs.
- Added verified reuse of complete runs with matching numerical fingerprints,
  NumPy, settings and observation balances, including compatible v0.5.1 outputs.
- Added equal-weight replicate means, sample SD/SEM, paired differences, final
  replicate tables, and Arial-first PNG/SVG ensemble plots. Partial ensembles are
  never silently averaged. All original per-run outputs remain available.
- Supplied the 30,000-a.u. MM32 ten-seed comparison, a small explicitly altered
  software diagnostic, tests and a step-by-step CLI walkthrough. Dependencies and
  the six numerical model modules remain unchanged.

## 0.5.1 - 19 October 2026 - Runtime reporting and typography

- Added elapsed wall time, current-process CPU time, initialisation/simulation/
  output phases, and actual solver throughput to run metadata and terminal output.
- Preserved timings on event limits and caught runtime failures; no clocks enter
  biological equations or random-number generation.
- Added one scoped Arial-first style for all static and animated outputs, font
  family/file controls, explicit fallback provenance, and strict font checks.
- Retained existing requirements and numerical model modules. Added a guarded
  module entry point and documented independent-run parallelisation as the next
  architecture step; this release does not launch workers.
- Existing v0.5.0 settings remain accepted; font defaults are supplied on replay.

## 0.5.0 - 19 October 2026 - Growth and signal production

- Added full mode: age-dependent triangular division appointments and per-cell
  BB170 Poisson signal production alongside binding, movement, and response.
- Added nearest-footprint daughter placement, retained parent attachments,
  newborn response clocks, cap/blocked-division handling, and explicit POP choices.
- Added CLI/JSON controls, lineage/population tables and plots, production balance,
  variable-population frame archives, and birth-aware GIFs with fixed colourbars.
- Added an MM32 full-mode configuration and a clearly labelled accelerated BB170
  growth demo; earlier diagnostic modes preserve their scope and event history.
- Includes v0.4.1's animation correction. Dependencies remain unchanged.
- Tests cover new rules and outputs. Published time-course comparisons and
  independent replicate averages remain assignment work to complete.

## 0.4.1 - 19 October 2026 - Stable animation geometry

- Fixed colourbar heights shifting during the first movement/response GIF frames.
- Settle automatic layout before recording and freeze the resolved axes geometry.
- Added a regression check for stable axes/colourbar bounds and numeric limits.
- Initialisation, kinetics, transport, expression, and saved numerical results are
  unchanged. This is a presentation fix, with no new biological interpretation.

## 0.4.0 - 19 October 2026 - Delayed expression and predicted light

- Added continuous per-cell expression under piecewise-constant receptor sensing,
  with analytical interval updates and no fixed-step approximation.
- Recorded the degradation-sign interpretation, branch crossing, initial expression,
  printed light transform, and population aggregation as RESP-01 through RESP-05.
- Added response mode, external rate/initial/midpoint/scale controls, and fixed
  light-map display limits. Numerical response is independent of sampling/logging.
- Added per-cell/population response tables, final/frame arrays, response curves,
  paired free-signal/light maps, and a second GIF. Existing movement views remain.
- Verified downstream response does not change stochastic transport/binding history.
- No growth, division, BB170 signal production, calibrated experimental light, or
  completed published-curve reproduction is claimed.

## 0.3.0 - 19 October 2026 - Binding and movement

- Added coupled free-signal diffusion, polymer carrier diffusion, and connected
  bacterial translations, with capacity/geometry checks and closed boundaries.
- Implemented the stable S22 mobility function and recorded obstruction,
  tethering, displacement, and D_B/m interpretations as MOVE-01 through MOVE-07.
- Added local propensity caches and Fenwick weighted event selection to handle
  frequent signal hops without fixed-step approximations.
- Added independent placement streams (INIT-05), explicit RNG scheme recording,
  and legacy placement replay for older files lacking the new setting.
- Added movement controls, a documented example, start/end maps, GIF animation,
  spatial arrays, position tables, compact/full logs, and per-event-type counts.
- Raised the default event cap to 500000; saved files retain explicit old caps.
- Added Pillow as an explicit animation dependency. No growth, signal production,
  expression, or luminescence is claimed.

## 0.2.0 - 19 October 2026 - Fixed-position binding

- Added all three reversible binding classes with Gillespie event selection and
  exponential waiting times, including individual polymer bridges and release.
- Implemented finite site capacities, signal conservation checks, and the S23
  polymer signal-site crowding cutoff.
- Added CLI/JSON binding controls, explicit polymer-rate requirements, independent
  dynamics RNG, trajectory/event tables, final arrays, and diagnostic plots.
- Preserved initial-only mode and its settings. Event-limited runs are explicitly
  incomplete and return exit code 3.
- Documented geometry, site allocation, fixed-anchor limitations, kinetic choices,
  and random-stream choices as BIND-01 through BIND-06. No movement, growth,
  signal production, expression, or luminescence is claimed.

## Documentation update - 19 October 2026

- Added the requested module, collaboration, and last-updated comment headers to
  all Python source and test files.
- Expanded comments explaining configuration precedence, grid indexing, placement,
  random allocation, output provenance, plotting, and validation checks.
- This update changes comments and documentation only; model rules, defaults,
  random draws, and executable behaviour are unchanged.

## 0.1.0 - Initial spatial state

- Added CLI and JSON settings; source edits are not required to change controls.
- Implemented seeded bacterial placement and initial free signal/polymer arrays.
- Added isolated output folders, position tables, plots, replay settings, source
  hashes, software versions, and run status metadata.
- Recorded paper sources, initialisation assumptions, and unresolved model details.
- Added initial-state validation and CLI integration checks.
- No stochastic trajectories or scientific reproduction claims yet.

