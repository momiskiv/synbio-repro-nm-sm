# Affinity controls — v0.7.0

The completed MM32 comparison showed lower P3 light in every matched seed pair,
while polymer-linked clusters developed. Affinity controls now test the model
consequences of enabling each binding interaction. They do not change the solver,
growth, bacterial signal receptors, expression rule or transport parameters.

## Source definitions and chosen kinetics

Supplementary page **S26**, “Modelled luminescence”, item 3, specifies P1 affinities
`K_PB=5, K_PS=0`, P2 `K_PB=0, K_PS=0.1`, and P3 `K_PB=5, K_PS=10`.
S21 defines affinity as forward rate divided by reverse rate. The individual
polymer rates are not supplied; our explicit reverse rates remain `0.001`.
Forward rates below implement those ratios without fitting them to a desired
light response. “P1-like” and “P2-like” identify affinity-based representations,
not a chemical reconstruction of the actual polymers.

| Condition label | Polymer entities | PB on / off | PS on / off | K_PB | K_PS |
| --- | ---: | --- | --- | ---: | ---: |
| `no_polymer` | 0 | 0.005 / 0.001 | 0.01 / 0.001 | — | — |
| `p1_like_1000` | 1000 | 0.005 / 0.001 | 0 / 0.001 | 5 | 0 |
| `p2_like_1000` | 1000 | 0 / 0.001 | 0.0001 / 0.001 | 0 | 0.1 |
| `p3_no_pb_1000` | 1000 | 0 / 0.001 | 0.01 / 0.001 | 0 | 10 |
| `p3_1000` | 1000 | 0.005 / 0.001 | 0.01 / 0.001 | 5 | 10 |

Unused polymer rates remain recorded in the no-polymer case so its request is
identical to the completed comparison. Zero affinity is implemented by a zero
forward rate with a positive reverse rate. All entities start unbound, so the
disabled class never acquires an attachment. Both bacterial signal binding and
polymer-dependent signal mobility remain active in every relevant condition.

**P2 is not P3 with only PB switched off:** its signal affinity is 100 times
weaker. The extra `p3_no_pb_1000` control preserves P3's signal affinity while
disabling bacterial binding. This keeps those two scientific questions separate.

## Start with the early-time pilot

From the v0.7.0 `polymer-clustering` folder, with your environment active:

```bash
python -m polymer_clustering --version
python -m polymer_clustering.batch --config examples/mm32_affinity_pilot.json --workers 2 --dry-run
python -m polymer_clustering.batch --config examples/mm32_affinity_pilot.json --workers 2 --output results/exercise5/mm32_affinity_pilot
```

Expect version `0.7.0` and **ten pilot runs**: five conditions, seeds 2026/2027.
The pilot retains the long-run lattice, 25 initial MM32 cells, 10,000 signals,
all transport/response/growth parameters and 1,000 entities in polymer conditions.
Only duration (1,000 a.u.), spatial frames (11), and event budget (10 million)
differ from the full plan. Sampling remains every 100 a.u. It is a shortened
scientific pilot, unlike the earlier tiny BB170 software diagnostic.

There should be no births before the default earliest division time of 5,000
a.u. This pilot examines early binding/sensing; it cannot establish the late-time
cluster/growth response. Two seeds provide an early consistency check, not the
ten-run average required for the full comparison. The 30,000-a.u. records cannot
be imported as 1,000-a.u. jobs; their requests/outputs represent different runs.

For the first review, send `batch.json`, `ensemble.png`, `replicate_summary.csv`,
`ensemble_trajectory.csv` and `pairwise_differences.csv` from the printed batch
folder. A small archive of those five files is sufficient for an initial review.
Keep the individual run outputs for detailed checks or animations.

## Read the controls and their contrasts

| Check or comparison | Interpretation |
| --- | --- |
| P1-like polymer-bound signals | Exactly zero at every observation: PS association is disabled |
| P2-like / P3-without-PB bonds and bridges | Exactly zero; polymer-connected largest cluster is one cell while cells are present |
| P3 minus P1-like | Enable PS binding at fixed PB affinity 5 |
| P3 minus P3-without-PB | Enable PB binding at fixed PS affinity 10 |
| P3-without-PB minus P2-like | Increase PS affinity from 0.1 to 10 with PB disabled |
| Each polymer condition minus no polymer | Compare the combined consequences of polymer presence, transport and its permitted binding |

In `pairwise_differences.csv`, values are **condition minus reference**, calculated
within each seed first, then averaged. Positive light differences mean the
condition has more light; negative values mean less. Each unordered condition
pair appears once, with the later-listed condition minus the earlier-listed one.
The original `paired_differences.csv` still compares every condition to the first
condition, no polymer. SD/SEM are descriptive across-seed statistics, and the
files do not run significance tests at every time point.

P1-like may still change signal transport despite having no signal attachments.
Disabling PB removes polymer bridges, but nearby daughters can still appear
adjacent on a map. The model's cluster metric uses explicit polymer connectivity,
not visual proximity. Enabling a binding class also changes subsequent spatial
histories; these contrasts measure the total model consequence of that control,
not an additive decomposition of individual reaction events.

The no-polymer control also removes polymer-dependent transport. It is not a
1,000-polymer condition with both affinities zero. A full factorial interaction
analysis would require that additional control. This first stage addresses
conditional PB/PS effects and the source P1/P2/P3 affinity definitions.

## Extend to the full matched comparison

`examples/mm32_affinity_full.json` specifies **50 planned runs**: five conditions,
ten seeds 2026–2035, duration 30,000 a.u., 61 frames, and the original event cap.
No-polymer and P3 requests are unchanged from the twenty runs already completed.

Set `PREVIOUS_BATCH` to the actual extracted/original **batch folder containing
`batch.json`**, not its tar archive or an individual seed folder. Replace this
example path and keep the quotes:

```bash
PREVIOUS_BATCH="/path/to/full_20261020T002045778556Z_521aaa9e"

python -m polymer_clustering.batch \
  --config examples/mm32_affinity_full.json \
  --workers 2 \
  --reuse-batch "$PREVIOUS_BATCH" \
  --output results/exercise5/mm32_affinity_full \
  --dry-run
```

Expect **20 reused and 30 pending**. Every imported run is checked against its
saved settings, numerical source, NumPy version, time grid and count balances.
Condition names do not substitute for those checks. Incompatible runs are not
silently reused; a source batch with no matching jobs is rejected. Duplicate
imports are rejected rather than counting a seed twice. Matching originals are
copied into the new batch and remain unchanged.

After reviewing the pilot and the full preview, use the same command without
`--dry-run` to start. The remaining thirty trajectories are substantial work;
their durations cannot be inferred reliably from the previous P3 runtime because
different affinities change event workloads. Use the already checked two-worker
setting if your memory and CPU allocation still permit it. Default remains one.

Resume an interrupted v0.7.0 batch with `--resume` and the printed batch folder,
as before. An old v0.6.0 batch must be resumed with its original package, or its
complete runs reused into a new v0.7.0 batch. Resume never changes an existing plan.

## Scope for the assignment

The primary question is whether PS binding contributes to reduced sensing/light,
and whether PB binding changes that response at matched PS affinity. Record the
outcome rather than requiring clustering to increase light. Source-matched
affinity ratios do not uniquely determine kinetics: fixing the off rate is still
a reconstruction choice. Polymer entity counts remain an uncalibrated proxy for
the publication's concentration axis, and outputs retain raw model light.

S26 also describes a 25-combination positive-affinity sweep using 0.1, 0.5, 1,
5 and 10 for both affinities, along with concentration variation. That larger
sweep and BB170 remain subsequent steps; this five-condition comparison is not
the full published surface or a claim of source-figure reproduction.
