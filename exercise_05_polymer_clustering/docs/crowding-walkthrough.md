# Matched PS crowding-sensitivity pilot

Version 0.9.0 tests whether PS-site occlusion contributes to the affinity result.
The intervention disables the PS cutoff while retaining P3 bacterial binding,
signal affinity and transport rules. This is a sensitivity experiment on the
reconstruction, not a new polymer chemistry or a published default. The default
remains **10 polymers per lattice site**.

## Design

| Setting | Standard P3 | P3 with PS cutoff disabled |
| --- | --- | --- |
| PS crowding threshold | 10 | `null` in JSON, `off` on CLI |
| PB forward/reverse rates | 0.005 / 0.001 | Identical |
| PS forward/reverse rates | 0.01 / 0.001 | Identical |
| Polymer count | 1000 | Identical |
| Strain, lattice, initial bacteria | MM32, 50 x 50, 25 | Identical |
| Seeds | 2026 and 2027 | Matched |
| Duration | 10000 a.u. | Identical |
| Sampling, saved frames | Every 100 a.u.; 21 frames | Identical |
| Growth and remaining controls | Existing source-default full-mode settings | Identical |

Four simulations are planned. The earlier 1000-a.u. pilot had no crowding at its
endpoints. This pilot extends to 10000 a.u., where earlier P3 time courses had
already begun to lose polymer-bound signals. Two seeds provide an exploratory
check rather than a ten-seed scientific ensemble.

Independent placement streams match initial bacteria, polymers and signals
within each seed. Event histories and realised clusters may diverge after the
rule is changed. The experiment holds rules and inputs fixed apart from PS
occlusion; it does not freeze the resulting geometry.

## Run

From the extracted v0.9.0 project folder, activate the existing environment.
Dependencies have not changed. First preview:

```bash
python -m polymer_clustering.batch \
  --config examples/mm32_crowding_pilot.json \
  --workers 2 \
  --dry-run
```

Expect four pending simulations. These are new 10000-a.u. requests. Earlier
30000-a.u. results remain separate and are not truncated or substituted.

```bash
python -m polymer_clustering.batch \
  --config examples/mm32_crowding_pilot.json \
  --workers 2 \
  --output results/exercise5/mm32_crowding_pilot
```

Two workers are a conservative starting point during concurrent writing on the
student's machine. This still requires simulation work: the few-second spatial
diagnostic runtime does not predict this batch's runtime. Record elapsed time;
a fixed completion-time estimate has not been benchmarked.

If interrupted, use the printed folder with the same package/environment:

```bash
python -m polymer_clustering.batch --resume "/printed/batch/folder" --workers 2
```

Review `ensemble.png`, `replicate_summary.csv`, `ensemble_trajectory.csv`,
`paired_differences.csv` and `batch.json`. The paired difference is **cutoff off
minus cutoff 10**. A negative light difference means stronger suppression after
removing PS occlusion. Population and signal-pool metrics help check alternative
interpretations. Bands show sample SD across two seeds, not confidence intervals.

## Endpoint density and capacity

Use the actual completed pilot folder:

```bash
CROWDING_BATCH="/printed/completed/batch/folder"

python -m polymer_clustering.spatial \
  --batch "$CROWDING_BATCH" \
  --conditions p3_cutoff_10 p3_no_ps_cutoff \
  --map-seed 2026 \
  --output results/exercise5/spatial_crowding_pilot
```

Send `spatial_summary.zip` and the batch summary files above. The ZIP supplies
endpoint evidence; the batch trajectory supplies the time course. Original
logs and animations remain in their run folders.

The diagnostic distinguishes the active cutoff from the fixed reference density
P >= 10. Its first comparison panel uses the fixed reference, so dense polymer
accumulation remains visible when PS occlusion is disabled. Maps label each
condition's active cutoff. Zero blocked slots when cutoff is off follows from
the intervention's definition; it alone is not evidence of a light effect.

## Interpretation and limits

If disabling PS occlusion increases polymer-bound signals and lowers receptor
occupancy and light, that supports a suppressive consequence of the intervention
in this model. Compare matched time courses and individual seeds before a longer
ensemble. Small or inconsistent effects would not support assigning the earlier
P3/no-PB difference entirely to PS occlusion.

The polymer-dependent signal-mobility function remains active with lambda=8.
Disabling PS occlusion does not disable all polymer crowding/transport effects.
Tethering, cluster movement, obstruction, PS release and response retain their
rules, though their realised trajectories can change downstream. This experiment
measures the combined downstream consequence of the PS rule change, without
separately identifying every localisation or transport mechanism.

## Parameter and provenance controls

`BindingConfig.ps_crowding_threshold` accepts a positive integer or `None`.
JSON uses `null`; CLI accepts `off` or `none`. Zero, negatives, booleans and
fractional thresholds are rejected. Explicit CLI `off` overrides a JSON
threshold; omitted controls retain JSON/default values.

`run.json` and `parameters.json` retain the actual setting. Altered rules are
marked CROWD-01 and labelled on individual/ensemble plots. Existing v0.8.0
configurations omit the field and retain cutoff 10.

Only `binding.py` changes among the six numerical modules. Default replay was
checked against v0.8.0 verification runs in fixed-position and coupled modes:
every CSV is byte-identical and all saved arrays match. Five numerical modules
and all pinned dependencies remain unchanged.

One explicit fingerprint bridge accepts completed legacy default-ten runs only
when the old/current binding hashes are the reviewed pair and the other five
numerical hashes and NumPy match. Missing legacy settings mean 10, never off.
Arbitrary edited sources are not accepted by this exception. Original settings
and files remain intact. Resume requires the original full package/environment;
cross-version reuse and read-only analysis are separate operations.
