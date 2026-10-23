# v0.11.1 interrupted-batch timing recovery

This small fix repairs the controller's final timing calculation after an
interrupted invocation. An invocation terminated by a reboot can lack
`wall_seconds`. The original controller directly indexed that missing field,
raising `KeyError` after the ensemble outputs had already been written.

Apply this fix only when all other batches using the same project have stopped
or completed. Keep the original v0.11.1 project and Python environment. There
are no dependency changes and no replacement simulation package to install.

## Apply and finish the standard batch

Extract the archive into the existing project directory. It creates the
`recovery-fix` subdirectory. From that project directory, set `BATCH` to the
existing standard 80-run batch folder, using the same path as the earlier
resume preview:

```bash
BATCH="/absolute/path/to/the/existing/standard/batch"

python recovery-fix/apply_resume_timing_fix.py --batch "$BATCH" --dry-run
python recovery-fix/apply_resume_timing_fix.py --batch "$BATCH"

python -m polymer_clustering.batch --resume "$BATCH" --workers 4 --dry-run
```

For the completed standard study, expect **80 retained; 0 pending**. Then finish
the batch record and regenerate its ensemble outputs:

```bash
python -m polymer_clustering.batch --resume "$BATCH" --workers 4 --progress on
```

With zero pending jobs, this performs verification and aggregation without
launching simulation trajectories. If the preview reports pending jobs,
confirm that you selected the completed standard batch rather than the next
study. Do not edit completion statuses manually.

## What changes

- The installer validates the exact reviewed v0.11.1 sources, original
  environment, scientific plan and completed run outputs before replacing files.
- It backs up the original `batch.py` and `batch.json` next to those files, using
  unique `.bak` names, and records the exact patch hashes and backup paths.
- Only the controller source fingerprint is updated in the batch record, so the
  normal strict resume checks can continue. No job, seed, scientific parameter,
  numerical module or individual completed-run file is changed.
- Unavailable controller invocation durations remain unknown. The full total
  `controller_wall_seconds_all_invocations` is JSON `null` when any invocation
  lacks a measurement. The separately named recorded-invocation total sums only
  measured durations, and the missing-invocation count makes the gap explicit.
- The prior `last_execution_wall_seconds` is preserved in the patch receipt and
  original backup. Execution-loop timing differs from the full controller time,
  which also includes verification, aggregation and rendering.
- Repeating the installer is safe. It can also complete a source-first update
  interrupted before the batch record was replaced.

The fix does not determine why the desktop returned to its login screen. Linux
logs are needed to distinguish a reboot from a desktop/session failure.

Fresh batches created after the fix use the patched fingerprint automatically.
If an older unfinished batch needs to be resumed later, run this same installer
for that batch while no other batch is active; it will update only its reviewed
controller identity after revalidation.

## Validation

Ten regression tests cover unknown/known/invalid timing, completed-result
preservation, environment/source/plan rejection, active-controller locking,
backup/idempotence, write-failure rollback, source-first interruption recovery,
and successful ensemble finalization with a missing old invocation duration.
Additional existing controller regression checks are recorded in VALIDATION.txt.
All six numerical module fingerprints are unchanged. Fixtures are software
checks and are not assignment data.
