# Model sources and reconstruction choices

Reference: **Lui**, L. T. et al. (2013). Bacteria clustering by polymers induces
the expression of quorum-sensing-controlled phenotypes. Nature Chemistry 5,
1058-1065. DOI: https://doi.org/10.1038/nchem.1793.

This project reconstructs published methods; it is not the authors' original code.
Version 0.5.0 implements initialisation, binding, transport, response, division, and signal production. Source pages below
refer to the numbered supplementary pages, not the main article pages.

## Implemented source mapping

| Feature | Source | Implementation |
| --- | --- | --- |
| 50 x 50 two-dimensional grid | S20, model detail 1 | `InitialConfig.grid_size` |
| 2 x 2 bacterial footprint; no bacterial overlap | S19-S20 | `initialise`, occupancy array |
| Signals and polymers can coexist with bacteria; multiple per site | S20, detail 1 | Separate integer count arrays |
| Initial bacteria: 25; model population cap: 80 | S23, detail 11 | Configuration defaults and validation |
| Initial signals: MM32 10000; BB170 0 | S21, detail 4 | `InitialConfig.signal_count` |
| Three reversible binding classes, mass-action propensities | S20-S21, details 2-4 | `reaction_channels`, `apply_reaction` |
| Gillespie direct event selection and exponential waits | S20-S21, detail 3 | `simulate_binding` |
| Polymer capacities: 2 bacterial ends, 5 signals | S21, detail 4 | Explicit per-polymer contact and signal arrays |
| Bacterial capacities: 120 polymer sites, 60 signal sites | S21, detail 4 | Equal per-site allocation: BIND-01 |
| Bacterium-signal rates: 0.001 forward, 0.002 reverse | S21, detail 4 | `BindingConfig` defaults |
| P1 affinities PB=5, PS=0; P2 PB=0, PS=0.1; P3 PB=5, PS=10 | S26, Modelled luminescence item 3 | `mm32_affinity_full.json`; fixed off-rate interpretation AFF-01 |
| Attached polymers can bridge neighbouring bacteria | S21, detail 4 | Individual contact records; geometry BIND-02 |
| Polymer signal-site cutoff at local P >= 10 | S23, detail 10 | `CROWDING_THRESHOLD`; interpretation BIND-04 |
| Cardinal movement; D_B=0.0005, D_P=0.01, D_S=0.05 per direction | S19-S22, detail 7 | `MovementConfig`, `MovementEngine` |
| Polymer-dependent signal mobility, lambda=8 | S22, detail 8 | Stable `signal_mobility` |
| Connected bacteria move as one unit with size-dependent rate | S19-S22, detail 9 | `clusters`, `cluster_can_move`, `move_cluster` |
| Moving clusters displace small chemicals | S20, S23 | Conservative chosen rule MOVE-04 |
| BSsense and continuous BSexpress response | S21-S22, detail 5 | `ResponseTracker`, `advance_expression` |
| Ksyn=0.05 and quoted degradation magnitude 0.01 | S22, detail 5 | `ResponseConfig`; sign choice RESP-01 |
| Light transform, Scale=2500, BSc=38 MM32 or 15 BB170 | S22, detail 6 | `light_from_expression`; RESP-04 |
| Age-dependent division, nearest empty daughter footprint, retained parent attachments | S19-S20 | `DivisionSchedule`, `nearest_empty_anchor` |
| Triangular lifetime mean 15000, width 20000; cap 80 | S23, detail 11 | `PopulationConfig`; POP-01 |
| BB170 signal production 0.005, random producer footprint site | S23, detail 12 | `produce_signal`; POP-05 |

Coordinates are lattice units, not micrometres. Arrays use `[y, x]` and the plot
origin is lower-left. Each bacterial coordinate is the lower-left occupied site;
its footprint occupies `(x, y)` through `(x + 1, y + 1)`.

## Explicit reconstruction choices

These choices are not claims about the original implementation.

| ID | Choice | Reason and limitation |
| --- | --- | --- |
| INIT-01 | Random sequential non-overlapping placement at integer anchors within closed edges | The precise initial placement algorithm is not provided in S19-S24. A shuffled list of all anchors is inspected until enough bacteria are placed. This is not uniform sampling over complete configurations. Dense placements can jam and fail explicitly. |
| INIT-02 | All signals and polymers start free, independently allocated uniformly across lattice sites | A documented initialisation assumption. This milestone does not assign receptor occupancy or gene expression. Initial binding states remain to be investigated. |
| INIT-03 | Explicit PCG64 generator with default seed 2026 | Our reproducibility choice, not an author-supplied seed. Same parameters and validated software reproduce numerical state arrays; cross-version identity is not promised. |
| INIT-04 | Polymer input is an entity count | The `polymers` option describes initial entities only. No conversion from physical concentration or the publication's `[P]` axes is established here. |
| INIT-05 | Default independent PCG64 placement streams SeedSequence([seed, 0, i]), with i=0 bacteria, 1 polymers, 2 signals | Added in v0.3 so changing polymer count preserves initial bacterial and signal positions. Some v0.2 seeded states intentionally change. `legacy_v02` preserves the former shared draw order; settings files lacking an RNG scheme are interpreted as legacy. |

Changing `grid_size` or `bacteria` from published defaults is a user override,
recorded in `parameters.json`; it is not labelled a matched published condition.
Using multiple seeds alone does not make these initial states completed simulation
replicates. Binding mode freezes positions. Movement mode couples reactions and
transport. Response mode adds expression and predicted light. Full mode also
adds growth and strain-specific signal production.

## Binding reconstruction choices

These describe fixed-position binding. MOVE-03 supersedes BIND-03 in movement
mode; the other binding assumptions remain active.

These choices define the current diagnostic system. Exact SSA samples that
specified system; it does not remove uncertainty about the source model.

| ID | Choice | Reason and limitation |
| --- | --- | --- |
| BIND-01 | Divide bacterial capacities equally across four occupied sites: 30 polymer sites and 15 signal sites per site | Preserves the stated 120/60 totals and local same-site mass action. The source does not specify this spatial allocation. Receptors do not redistribute between sites in this version. |
| BIND-02 | First polymer attachment requires a bacterium at its anchor. A second end can bind a site on another bacterium sharing an edge with any site of its current host. Duplicate target sites are counted once; diagonal contact alone is excluded. One polymer connects at most two distinct bacteria. | Explicit interpretation of 'neighbouring bacteria'. Per contact propensity is forward rate times free polymer ends times free bacterial sites. Same-bacterium double attachment is excluded; this affects kinetics and requires later sensitivity analysis. |
| BIND-03 | Freeze all coordinates, including each polymer's original anchor, throughout the binding diagnostic | Deliberately staged subsystem, not the full paper model. Polymer-signal reactions and release stay at the anchor even if a polymer remains attached only to a neighbouring bacterium after its first contact releases. Subsequent bridging uses the remaining host's neighbours. No molecular geometry, polymer relocation, or cluster transport is inferred. These choices must be revisited when movement is implemented. |
| BIND-04 | Use all polymer entities at the anchor for the local P >= 10 signal-site cutoff. Suppress new PS binding while preserving existing PS attachments and their ordinary release. Recompute accessibility from current counts. | S23 gives the free-signal-site cutoff. It does not fully specify counting of attached polymers or reversibility of occlusion. Polymer counts cannot change locally in this fixed-anchor stage. S20 mentions loss of cell-binding sites too, but supplies no corresponding quantitative PB rule; no extra PB cutoff is invented. |
| BIND-05 | Require explicit polymer forward and reverse rates. Example: PB 0.005/0.001, PS 0.01/0.001, hence ratios 5 and 10 | S26 reports P3 ratios 5 and 10, but not these individual rates. These are chosen kinetics, not estimated or fitted values. The example's 100 entities and 100 a.u. duration are demonstration settings, not a published condition. BS defaults are the published 0.001/0.002 pair; overrides are recorded as user-selected kinetics. |
| BIND-06 | Give dynamics its own PCG64 stream via SeedSequence([seed, 1]); keep the original placement stream unchanged | Reproducibility choice. Observation times are independent of events, and the final requested time is included. A maximum event count stops with explicit incomplete status. Full channels are rebuilt and conservation/capacity checks run after each event; performance optimisation remains future work. |

Signal accounting is `free + bacterium-bound + polymer-bound = initial total`.
Each PB bond occupies one bacterial site and one end of exactly one polymer.
The unattached-plus-attached polymer count is also conserved; a bridge still
represents one polymer. Signal binding and bacterial binding use separate sites
on a polymer. PS sites remain available on attached polymers unless occluded by
the stated crowding rule. No receptor-bound signal is counted as free.

## Movement reconstruction choices

| ID | Choice | Reason and limitation |
| --- | --- | --- |
| MOVE-01 | Four orthogonal directions with closed edges. A free entity at a site contributes D per available direction; aggregated signal hop propensity is free count times D_S times f(P) times the edge multiplier. | S22 gives propensity to a specific neighbour. Grouping indistinguishable signal hops preserves these rates. No division by four and no periodic wraparound are introduced. |
| MOVE-02 | An edge touching a bacterial footprint at either endpoint receives `obstruction_factor`, default 0.1, for signals and independently moving polymers. | S19 describes slower obstructed directions, while S22 also calls some directions unavailable. The exact lattice geometry and reduction are unspecified. This symmetric, permeable interpretation permits signals to enter and leave receptor-containing sites. 0 gives strict blocking; 1 removes this obstruction. These alternatives require later comparison. |
| MOVE-03 | A polymer without bacterial contacts diffuses at D_P and carries its bound signals. A tethered polymer has no independent hops; its anchor is the first surviving contact. Releasing the last end leaves it at its last anchor. BS-bound signals move only with their bacterium. | Explicit carrier/coordinate convention, replacing the frozen-anchor diagnostic choice BIND-03. Contact reassignment can relocate a polymer and its bound signals within the adjacent-cell geometry; it is part of the PB event. No bound signal is also counted in free diffusion. Polymer-carried signals do not alter D_P in this version. |
| MOVE-04 | For a one-site cluster translation, shift free signals and untethered polymers originally in the union of old/new bacterial footprints one site in the same direction. Bound signals and tethered polymers translate with their owners. Pushed free entities at an outer edge remain at that edge. Counts at shared destinations add. | Concrete, conservative interpretation of S20/S23 displacement, not a recovered original algorithm or a hydrodynamic calculation. Updates are simultaneous to avoid double pushing. Remote free chemicals are unchanged. Free chemicals may share bacterial sites. Anchored polymers outside the moving group are not pushed. The `advection` control supports later assessment against no displacement. |
| MOVE-05 | Use D_B/m per direction for an m-bacterium connected component, including m=1. Reject a translation if any footprint hits an edge or another component. Translate all members and contacts together; rebuild components after PB events. | S22 prints D_S/m but describes its numerator as the bacterial rate. Using the separately specified D_B=0.0005 is a documented notation correction. Connectivity comes from explicit polymer bridges, not mere bacterial adjacency. |
| MOVE-06 | Evaluate f(P) = (tanh((lambda-P)/2)+1.2)/2.2 with lambda=8 and P equal to the current local count of all anchored/mobile polymers. | Algebraically matches the S22 exponential expression without overflow. It is not renormalised: f(0) is approximately 0.999695, and its high-P limit is 1/11. Local counts are a concentration proxy under INIT-04, not a physical concentration calibration. Use source-site P for each outgoing hop, as the stated local rule; spatially varying transport may therefore be asymmetric across an edge. |
| MOVE-07 | Hierarchical exact SSA with cached site/polymer rates and Fenwick weighted indices; local checks per affected site and full checks at observations, every `validate_every` events, and at the end. Rebuild geometry after PB events or cluster moves, and refresh partial sums periodically. Dynamics uses PCG64(SeedSequence([seed, 2])). | Software organisation and performance improvement, with no fixed-step approximation. Cache rates are checked against full reconstruction in tests. Compact logs omit signal-hop rows but keep total per-type counts; full logs remain selectable. Sampling, frame count, and logging choices do not change event draws. Changing validation interval can change floating-point summation history, so replay uses all saved settings. |

Movement/response diagnostic modes conserve the initial total of signals and
polymer entities. Full mode conserves initial plus produced signals, and all
polymer entities. Cluster movement preserves all bacterial
footprints, bound signal pools, and PB endpoints. Signal-site crowding is
recomputed after polymer relocation, so accessibility can recover when local
polymer count falls below 10; existing PS bonds survive the cutoff and can release.

## Response reconstruction choices

| ID | Choice | Reason and limitation |
| --- | --- | --- |
| RESP-01 | Use Ksyn=0.05 and positive Kdeg=0.01 with dE/dt=-Kdeg*E in the falling branch | S22 prints a leading minus sign and Kdeg=-0.01, which together produce growth. Taking the quoted magnitude with the printed negative derivative preserves the described decay. This is an explicit interpretation, not a verified author-code correction. User overrides are saved. |
| RESP-02 | Integrate exactly over intervals of constant per-cell receptor occupancy S. Rise: E(t+dt)=S+(E(t)-S)*exp(-Ksyn*dt). Fall: max(S,E(t)*exp(-Kdeg*dt)); equality is stationary. | Binding/unbinding changes S discretely. In the falling branch, reaching S switches to zero derivative; decay must not continue below the unchanged sensing level. This avoids branch overshoot caused by a coarse update. The original implementation's integration scheme is unknown. There is no feedback from expression/light into SSA propensities, so the event clock remains unchanged. |
| RESP-03 | Initially E=0 for every cell by default; permit a recorded scalar override between 0 and 60 | Source initial expression is not specified. Initial sensing is determined by INIT-02, so it starts at zero. E is a continuous proxy in receptor-occupancy units, not a simulated integer protein count. Cells retain their expression state during movement. Full-mode daughter initialisation is specified separately in POP-04. |
| RESP-04 | Retain lum=Scale*(1-(exp(delta)-1)/(exp(delta)+1)), delta=0.25*(BSc-E), implemented stably as 2*Scale*exp(-logaddexp(0,delta)) | This preserves the printed factor of two, the midpoint lum=Scale, and a nonzero baseline at E=0. There is no subtraction, clipping to Scale, or imposed zero threshold. The theoretical asymptote is 2*Scale, although E is bounded by 60 here. BSc defaults to 38 for MM32 and 15 for BB170; Scale defaults to 2500. Arbitrary model light is not calibrated to photons or instrument units. |
| RESP-05 | Compute light per bacterium, then save raw sums and means; no fold-change or population normalisation is inferred | The nonlinear transform of mean expression would not equal mean per-cell light. The source's final population normalisation remains insufficiently specified. Empty-population sums and means are set to zero. Display maps repeat each per-cell value over its footprint but aggregation counts each cell once. Automatic per-run colour scales are labelled; `light_map_max` permits fixed visual comparisons without changing numerical outputs. |

Response calculations use per-cell clocks updated only at receptor-occupancy
changes. Observations evaluate the closed-form solution without modifying these
clocks. Adding frames or sample points cannot change expression integration.
Other events, including cluster movement, preserve each bacterium's sensing total
and expression identity. Each response-enabled simulation therefore has the same
binding/movement history as its response-disabled counterpart under identical
initial, kinetic, transport, and random settings. These are conditional model
predictions, not a completed reproduction of experimental luminescence curves.

## Population reconstruction choices

| ID | Choice | Reason and limitation |
| --- | --- | --- |
| POP-01 | Symmetric triangular cell lifetimes with mode/mean 15000 a.u. and full width 20000: support [5000,25000]. All initial cells start age zero. Reset the parent clock and start a daughter clock after successful division. | S23 gives mean and width but not endpoints, symmetry, or initial ages. Lifetimes are sampled once per cycle, not implemented as a constant exponential division propensity. Width zero is an exposed deterministic diagnostic. Initial-age sensitivity remains unresolved. |
| POP-02 | Search all empty in-bounds 2 x 2 footprints. Minimise squared Euclidean anchor distance to the current parent position; choose uniformly among tied minima. | Concrete interpretation of nearest empty space in S19. Bacteria cannot overlap; pre-existing free chemicals may share the daughter's sites. Parent coordinates, all bound pools, polymer contacts, and expression are retained; no bound material is partitioned. |
| POP-03 | Cap population at 80, with an exposed lower cap for diagnostics. If no valid footprint exists, retry after 100 a.u. without drawing a fresh lifetime until a successful birth. | Source specifies the cap but not blocked-division behaviour. No parent moves, material disappears, or daughter is placed remotely outside the domain. Each blocked attempt counts toward the event limit and is logged. |
| POP-04 | Daughter has zero bound molecules and default expression zero, independently configurable from initial-cell expression. Its response clock starts at birth. | Daughter attachments are specified by S19; expression inheritance is not. No dilution or partitioning of parental expression is introduced. This choice needs biological sensitivity assessment. |
| POP-05 | BB170 produces at 0.005 per cell per a.u. as a homogeneous Poisson process between population changes; MM32 produces zero. A producer is selected uniformly, then one of its four sites uniformly. | Interpretation of S23's fixed production rate as a per-cell stochastic propensity. Overrides for either strain are recorded experiments. No signal degradation is specified. Every production event adds one to both the free pool and the conserved total target. |
| POP-06 | Scheduled division interrupts the exact chemistry/transport/production SSA. Resample exponential waits after a division; fixed-rate hazards are memoryless. Lifetimes and placements use independent PCG64 streams [seed,3] and [seed,4]. | Solver/reproducibility choice. Equal-time appointments use cell-ID order, precede stochastic events, and are all processed before a same-time observation, unless the event cap interrupts them. Rebuild geometry after a successful birth; production refreshes the affected site's channels. Samples and frames never redraw division lifetimes. |
| POP-07 | Permit mean equal to half-width, with positive mean: S42 mean 10000 and retained S23 width 20000 give support [0,20000]. Resample finite-precision zero draws so scheduled lifetimes stay positive. | A continuous triangle has no atom at zero. S42 does not restate the width; retaining 20000 is an explicit fixed-width interpretation. Negative support remains invalid. Old positive-support requests keep exactly the previous values and RNG consumption. Implemented in v0.11.1 for BB170 collection. |

Full mode includes these processes without establishing a unique reconstruction of
the original implementation. The accelerated BB170 example changes mean lifetime
to 40 and width to 20 a.u. for demonstration, explicitly distinct from the source.
Default MM32 initial signals remain 10000; BB170 starts at zero and accumulates
production. The default cap applies to both strains, and no cell death is added.

## PS crowding sensitivity (v0.9.0)

**CROWD-01:** `ps_crowding_threshold` defaults to 10, retaining the S23
interpretation. An alternate positive threshold or null/off is a labelled
sensitivity intervention, applied only to new PS on channels. Already-bound
signals retain ordinary off channels. PB sites/rates and the separate lambda=8
signal-mobility function are unchanged. The four-run pilot retains PB and changes
this single rule with matched seeds and all other controls identical. Downstream
event histories and geometry may change; no new polymer chemistry is implied.

The known fixed-ten binding fingerprint is bridged only to the reviewed
parameterised source and only for default-ten requests. All other numerical
fingerprints and NumPy must match. Actual default verification replay matched
every saved CSV and array. The exception does not accept arbitrary edited sources
or import old default-ten runs as cutoff-off observations.

## Endpoint spatial diagnostics (v0.8.0)

The separate `polymer_clustering.spatial` command reads final arrays without
running the solver or altering original files. Its measurements use the existing
P >= 10 PS cutoff and five PS slots per polymer. Polymer localisation means the
anchor's overlap with a bacterial footprint, not proof of a current tether.
Empty PS slots at crowded anchors are counted separately from retained bound
signals. "Without cutoff" propensities are hypothetical rates at the same saved
state, not altered-rule simulations. Existing frame arrays lack polymer
coordinates over time; endpoint associations cannot isolate a causal mechanism.
See spatial-walkthrough.md for all definitions and limitations.

## Affinity-control reconstruction choices (v0.7.0)

| ID | Choice | Reason and limitation |
| --- | --- | --- |
| AFF-01 | Keep both polymer reverse rates at the existing chosen 0.001 and obtain affinities by changing forward rates. P1: PB 0.005, PS 0; P2: PB 0, PS 0.0001; P3: PB 0.005, PS 0.01. | Matches S26's stated ratios but does not recover unpublished individual kinetics. Holding reverse rates fixed defines one reproducible sensitivity experiment; other equal-ratio kinetics can change time courses. |
| AFF-02 | Include a separate P3-without-PB control (PB forward 0, PS forward 0.01). | S26's P2 has K_PS=0.1, not 10. This additional chosen ablation changes PB alone relative to P3, rather than conflating loss of bacterial binding with weaker signal affinity. |
| AFF-03 | Use the existing 1,000-entity MM32 setting, identical remaining controls and matched seeds. Start with two seeds to 1,000 a.u.; full comparison uses ten to 30,000 a.u. | Staged investigation of our reconstruction, not the full published affinity/concentration sweep. The pilot stops before any default division. No new calibration of polymer count, luminescence normalisation, geometry or kinetics is implied. |

Zero-affinity forward reactions stay disabled from the unbound initial state.
Polymer-dependent mobility still uses the local entity count for every polymer
condition. P1-like is therefore not a pure clustering-only perturbation, and
P2-like is not a removal of all polymer-associated transport effects. Conditional
P3-minus-control contrasts measure the total consequence of a changed binding
class, including downstream spatial/transport changes. A both-affinities-zero
polymer condition would be required for a full two-by-two binding interaction
analysis at fixed polymer presence. No such interaction claim is made here.

These AFF labels describe the experiment configurations and are documented here;
the existing single-run solver's interpretation list remains unchanged. Numerical
source modules and reaction rules are unmodified.

## Concentration collection - version 0.10.0

| ID | Choice | Rationale and limits |
| --- | --- | --- |
| CONC-01 | Sample MM32 P3 at 0, 1000, 3000 and 7000 initial entities, seeds 2026/2027, to 30000 a.u. every 100; retain cutoff ten and existing full affinity settings, including 61 frames. | S26 explicitly states 0-7000 in steps of 1000 and time 0-30000. Some legends/captions conflict (70000 and 0-3000); follow the explicit range and record the ambiguity. This focused two-seed pilot is not the full sweep or its ten-run average. Four compatible complete control/P3-1000 runs can be imported. Counts are not calibrated experimental concentrations. |

SI S31-S39 motivate comparing time-course shape and fold change across loadings.
Suppression and enhancement are possible outcomes, not targets to fit. No
biological equation, rate, geometry, random stream, growth rule, initial
expression or normalisation changes in v0.10.0. Binding/movement files gain
read-only progress hooks and reviewed fingerprint compatibility exceptions.

Before later BB170 collection, evaluate RESP-03's initial-expression assumption.
With non-negative expression and E(0)=0, the monotonic S22 light transform starts
at its minimum. It cannot display the initial decline in the source's BB170
example curves. E(0) is unspecified in the source; do not silently change the
assumption or commit a long ensemble before diagnosing the comparison.

## Unresolved issues for the complete model

BB170 Part 3 plans are documented in [the collection walkthrough](bb170-walkthrough.md).
They cover S40-S42's standard, low-density and faster-growth regimes at the
explicit S26 loading range, ten seeds each. Initial and daughter expression
remain zero; collection does not resolve the following scientific ambiguities.

- S21 gives some binding affinities as ratios without fully specifying both
  forward and reverse kinetic rates. Equal ratios need not imply equal dynamics.
  BIND-05 exposes this uncertainty rather than resolving it scientifically.
- S21-S22 define signal sensing and expression, but their initial values need
  an initialisation choice. RESP-03 records E=0; its effect needs evaluation.
- S22 prints a negative degradation parameter and an equation with an additional
  minus sign. Literal substitution produces growth rather than the described decay.
  RESP-01 records the decay interpretation; the authors' code is not available here.
- S22 prints `D_S/m` for cluster movement while describing the numerator as the
  bacterial rate, elsewhere named `D_B`. MOVE-05 records the chosen correction;
  it is not confirmation of the authors' implementation.
- S19 describes reduced movement through obstructed directions; S22 also discusses
  unavailable directions. MOVE-02 is an explicit provisional rule, requiring
  sensitivity analysis rather than a claim of unique reconstruction.
- S23's crowding cutoff disables free signal-binding sites locally at `P >= 10`;
  BIND-04 defines the diagnostic choice. Movement updates accessibility using
  current local populations; this reversibility assumption needs evaluation.
- S23 gives triangular division times with mean 15000 and width 20000. Symmetry,
  endpoints, starting ages, and scheduling are recorded in POP-01/POP-06; the
  original implementation and consequences of these choices remain unverified.
- S26 states time 0-30000; multiple captions say 0-3000. Main Figure 2 labels a
  polymer range up to 70000, whereas S26 specifies 0-7000. Check targets individually.
- BIND-01 defines a provisional spatial site allocation and MOVE-04 a provisional
  cluster transport rule. Their effects require evaluation. RESP-03/RESP-05 define
  initial expression and raw aggregation; the publication's normalisation remains
  unresolved. No claim of matching experimental amplitude is made.

## Implementation improvements versus biological changes

CLI/configuration controls, typed functions, output isolation, metadata, version
recording, source hashes, and validation improve software reproducibility. They
do not demonstrate improvement of the scientific model. Any later change to
biological rules requires a separate entry, rationale, and comparison.

