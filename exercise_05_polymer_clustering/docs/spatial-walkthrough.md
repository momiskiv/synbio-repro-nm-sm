# Endpoint polymer distribution and crowding

In v0.9.0 the diagnostic uses each run's saved PS cutoff (missing legacy value
means 10; null means off). It also records `sites_at_source_density_threshold`
and `polymer_dense_fraction_at_10`, both evaluated at the fixed density P >= 10
independently of whether occlusion is enabled. The first plot uses this fixed
reference density. Maps label the active cutoff on each row. The original
v0.8.0 definitions below describe runs using the default threshold ten.

Version 0.8.0 adds **post-processing of completed runs**. It reads the saved
final polymer identities, anchors, contacts and signal counts. No SSA event,
new seed, worker process or simulation is executed. The six numerical model
modules and all dependencies are unchanged.

This closes an observational gap in the affinity comparison: whether the final
P3 and P3-without-PB states differ in polymer distribution, PS-site accessibility,
and signal capture. It does not by itself identify the causal contribution of
crowding to the entire light trajectory.

## Run against your existing affinity batch

From the extracted v0.8.0 project directory, activate the existing environment.
Set the variable below to the **completed batch folder containing `batch.json`**.
Use the original results directly; do not copy the 1-GB collection into the new
project or rerun the simulations. The placeholder needs your actual path.

```bash
AFFINITY_BATCH="/path/to/completed/affinity/batch"

python -m polymer_clustering.spatial \
  --batch "$AFFINITY_BATCH" \
  --conditions p3_no_pb_1000 p3_1000 \
  --map-seed 2026 \
  --dry-run
```

The preview should validate **20 completed runs** at t=30000 a.u. If a source is
missing, incomplete, changed, incompatible or outside the batch, the command
stops without creating output. The loader verifies recorded settings, numerical
source/NumPy compatibility, observation balances, final array geometry/capacity,
and agreement of final arrays with final counts. Compatible v0.5.1/v0.6.0/v0.7.0
completed runs are accepted. The package version alone is not used as a proxy
for numerical compatibility.

Then run the same analysis without `--dry-run`:

```bash
python -m polymer_clustering.spatial \
  --batch "$AFFINITY_BATCH" \
  --conditions p3_no_pb_1000 p3_1000 \
  --map-seed 2026 \
  --output results/exercise5/spatial_affinity
```

No `--workers` option is needed. The terminal prints the analysis duration,
unique output folder and a **`spatial_summary.zip`** to send for review. This is
much smaller than the simulation archive. It contains CSV summaries, provenance
and two PNG figures, without event logs, animations or the full site table.

Omit `--conditions` to include all five conditions. Add `--seeds 2026 2027` to
inspect a smaller matched subset. All selected conditions must contain the same
seed set. `--map-seed` only selects the illustrative map, not the ensemble seeds.
Arial remains the default; the usual `--font-path`, `--font-family` and
`--strict-font` controls apply and the resolved font is recorded.

## Read the outputs

| File | Interpretation |
| --- | --- |
| `spatial_replicates.csv` | One endpoint row per condition/seed; retained sensing/light counts plus new spatial metrics |
| `spatial_ensemble.csv` | Equal-weight run means, sample SD and SEM; not pooled cells or a confidence interval |
| `spatial_paired_differences.csv` | Endpoint differences within each seed, then summarised; later-listed condition minus earlier-listed condition |
| `spatial_sites.csv` | Full auditable site table for all selected runs, with coordinates, density, occupancy and PS accessibility |
| `spatial_comparison.png/svg` | Four endpoint comparisons; each dot is one seed, bars show mean ± sample SD |
| `spatial_maps.png/svg` | One selected seed across conditions, with identical geometry and shared colour limits per column |
| `spatial_report.json` | Analysis scope, input hashes, original settings/versions, font and elapsed/CPU time |
| `spatial_summary.zip` | Compact upload bundle; omits the full site table and SVGs, which remain available locally |

For the suggested condition order, paired contrasts are **P3 minus P3-without-PB**.
A positive light difference means PB increases light at the endpoint; a positive
blocked-slot difference means more unoccupied PS slots are inaccessible in P3.
The smallest selected seed is mapped by default. It is an example, not a
representative or averaged spatial configuration.

## Metric definitions

PB uses polymer ends and PS uses separate signal slots. Bacterial attachment
does not directly consume PS slots. It can change where polymers accumulate,
how they move, local signal mobility and the availability of empty PS slots.

| Metric | Definition and units |
| --- | --- |
| `max_local_polymers` | Largest count of polymer anchors at any lattice site |
| `polymer_occupied_sites` | Number of sites containing at least one polymer |
| `crowded_sites` | Number of sites with local P >= 10, the existing PS cutoff |
| `polymers_at_crowded_sites` | Polymer entities at those sites; distinct from the number of sites |
| `polymer_crowded_fraction` | Fraction of all polymers located at crowded sites |
| `polymers_on_bacterial_sites` | Polymer anchors on bacterial footprints, whether or not currently tethered |
| `polymer_on_bacteria_fraction` | The preceding count divided by the polymer count; not the attached fraction |
| `empty_ps_slots` | Sum of five minus bound signals for each polymer |
| `blocked_empty_ps_slots` | Empty PS slots on polymers at crowded sites; previously bound signals are retained |
| `blocked_empty_ps_fraction` | Blocked empty slots divided by all empty PS slots |
| `ps_capacity_occupancy` | Polymer-bound signal count divided by five times the polymer count |
| `ps_bound_on_bacterial_sites` | Polymer-bound signals carried at anchors on bacterial footprints |
| `ps_bind_propensity` | Current aggregate PS on-event rate, in events per model a.u., using the existing crowding cutoff |
| `ps_bind_propensity_without_cutoff` | Hypothetical aggregate PS on-event rate with all empty slots accessible at this *same fixed state* |
| `ps_bind_propensity_lost_to_cutoff` | Difference between those two instantaneous rates, not lost binding events over time |
| `ps_unbind_propensity` | Current off rate times total polymer-bound signals, in events per model a.u. |
| `mean_free_signals_bacterial_sites` | Free signal count averaged over bacterial footprint sites |
| `mean_free_signals_other_sites` | Free signal count averaged over all other lattice sites |

Fractions are stored from 0 to 1 in CSVs and displayed as percentages in the
plots. A zero denominator produces an undefined blank, not a measured zero.
The corresponding ensemble has n=0 if every value is undefined. One replicate
has a mean but undefined SD/SEM. Actual zero counts remain valid zeros.

## What would fill the explanatory gap?

If P3 shows higher local density, more blocked empty PS slots and reduced PS
capacity occupancy across matched seeds, that would support a spatial crowding
explanation for its weaker sequestration. The maps show where those differences
occur. If the crowding difference is small, inspect polymer localisation and
local free-signal distributions as other explanations.

These observations do not separate tethering, clustering, mobility and crowding
as causal mechanisms. A future intervention could change one rule at a time,
with matched seeds and explicitly labelled reconstruction choices. The offline
"without cutoff" propensity is not that intervention: all coordinates and
signal counts are held fixed during this calculation.

The existing `movement_frames.npz` retains free signals and bacterial positions,
and `response_frames.npz` retains cellular response. Neither contains a polymer
coordinate time series. No polymer/crowding time course or animation is inferred
from the final arrays. No physical concentration calibration is introduced.

Original simulations remain unchanged. Use the original package/environment to
resume an interrupted simulation batch; this new package can analyse compatible
completed runs without changing their original software identity.
