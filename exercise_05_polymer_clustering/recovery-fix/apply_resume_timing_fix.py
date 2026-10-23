# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Apply the reviewed controller-only v0.11.1 interruption recovery fix.

Run from the original project directory with its original Python environment:
    python recovery-fix/apply_resume_timing_fix.py --batch /path/to/batch --dry-run
    python recovery-fix/apply_resume_timing_fix.py --batch /path/to/batch

The installer validates completed runs, preserves backups, and records the exact
controller source change. It never changes jobs, seeds, scientific parameters,
numerical modules, completed run files, or an unknown elapsed time into zero.
"""

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import stat
import sys
import tempfile
import uuid


PATCH_ID = "v0.11.1-resume-timing-20261022"
# Filled from the reviewed v0.11.1 release and the tested replacement below.
RELEASE_SOURCE_SHA256 = {'__init__.py': 'e1de1b2c0378b68307738464926785e7eb8ad19403d2fc967a9740980dd122b2',
 '__main__.py': 'ca7722b7f890497e8c12333406f540c2b3337b0d6aa141e33fe348d0b9b0161e',
 'batch.py': '8c372950cf1c89690bd7677ad35f83f15fa806987f28ae71c9b9abf46a46312f',
 'batch_output.py': '0147ac6aa17462f549b444a9c192cd308f0d256095672e0fe4cfe4a0dced5d5e',
 'binding.py': '1c5f02f0af61df5f05c10314fdcea5c7f2b767f7d2fca10b45c63040dcc7443f',
 'binding_output.py': '74064172c7d708cbc0db99c46ca7af6af6e8ea4159bac871b22e7ce59fd6c351',
 'cli.py': 'f3c7aa6b56d6a0048cb50714ba5494cd2ea3e73b99adaa3b013075423f3e6206',
 'model.py': 'eb2f62350d7ca61e30cb371899e4b3968706e5161ce1445dcf2b4907ff186339',
 'movement.py': '7c49b22443bb309b1efbd8fdd0053b2e7f2190e27854706f770b5d7bfe68ea73',
 'movement_output.py': 'dd70aa143dbb5b733c41cf915a265449f5a5fb14adbd93b1a17d94ce7dd43468',
 'output.py': '8c4f4a666b7884583c288c7e2ec92c68d7f21d5e854f59a49add2cec4fcb0c3e',
 'plotting.py': '5ea3ae4b6c97f95c3f04fa7a589717dd5792eea87765a5b5a629e76dddced82e',
 'population.py': 'cbb707ec8f375ed4382cb7402d09e9e89e3cdd2cb57e94efb0f818407df9422a',
 'population_output.py': 'fa042cf4ad86663e504064ad2a5d2189f171883c46b89835b80ae867fe1b0896',
 'progress.py': '9100e5fec7ea57e48e325b819f41119a934b6edd9c57344515ce8c25eaf0a1e8',
 'rates.py': 'ccd3f274c02441c51729a069af645019f58d9c8a5457ba0a4a67e242e6c146b5',
 'response.py': 'c19faac0756296207a36c317c3a4a54f0b8e82d49ba471cb9718175225b1d115',
 'response_output.py': '72d27c2f36a8bdf4341d24cbfe1191c0a2b08f5f696f5157cb872ba7bf4ad62b',
 'spatial.py': 'f525aa07999fda6c73bd0c317f732e0c179fd89942fa5eecce542ccd7d7d06bc',
 'spatial_output.py': '2e3399e8aa3771a827e8aa6d873159611c815eaa5e527f853e7a1eaa7194968b',
 'timing.py': 'd45f5a1bfcaf0e4f727fdbdaaea2f475d10d04fd7dfce09c3891d774d6274fe3'}
PATCHED_BATCH_SHA256 = "356a9a3fb1b814f3050138cefa3b81fba5e4041fd9a83f8fba736cdea5c0a65b"

OLD_TIMING_BLOCK = '''        manifest["controller_wall_seconds_all_invocations"] = sum(
            i["wall_seconds"] for i in manifest["invocations"]
        )
'''
NEW_TIMING_BLOCK = '''        record_controller_wall_time(manifest)
'''
HELPER_ANCHOR = "def finish_batch(directory: Path, manifest: dict, imports: list, workers: int,\n"
TIMING_HELPER = '''def record_controller_wall_time(manifest: dict) -> None:
    """Keep interrupted invocation times unknown while summing measured times.

    A reboot can prevent the finally block from recording an elapsed duration.
    UTC start/finish labels cannot recover a monotonic-clock measurement. Keep
    that gap explicit instead of inventing zero seconds or including downtime.
    """
    recorded = []
    missing = 0
    for invocation in manifest["invocations"]:
        elapsed = invocation.get("wall_seconds")
        if elapsed is None:
            missing += 1
            invocation["timing_status"] = "unavailable_after_interruption"
        elif (type(elapsed) not in (int, float)
              or not math.isfinite(elapsed) or elapsed < 0):
            raise ValueError("Invalid recorded invocation wall time.")
        else:
            recorded.append(elapsed)
    measured = math.fsum(recorded)
    manifest["controller_wall_seconds_recorded_invocations"] = measured
    manifest["controller_wall_seconds_missing_invocations"] = missing
    manifest["controller_wall_seconds_all_invocations"] = None if missing else measured


'''


def digest(data):
    """Fingerprint exact bytes, including comments, for a reversible patch."""
    return hashlib.sha256(data).hexdigest()


def replacement_source(original):
    """Require the exact reviewed source and both unique replacement anchors."""
    if digest(original) != RELEASE_SOURCE_SHA256["batch.py"]:
        raise ValueError("Controller is not the reviewed unmodified v0.11.1 file.")
    text = original.decode("utf-8")
    if text.count(OLD_TIMING_BLOCK) != 1 or text.count(HELPER_ANCHOR) != 1:
        raise ValueError("Expected timing code was not found exactly once.")
    updated = text.replace(OLD_TIMING_BLOCK, NEW_TIMING_BLOCK).replace(
        HELPER_ANCHOR, TIMING_HELPER + HELPER_ANCHOR
    ).encode("utf-8")
    if digest(updated) != PATCHED_BATCH_SHA256:
        raise ValueError("Replacement controller fingerprint does not match the reviewed fix.")
    compile(updated, "batch.py", "exec")
    return updated


def verify_software(previous, current):
    """Permit only this exact controller patch; retain all environment checks."""
    for software in (previous, current):
        expected = dict(RELEASE_SOURCE_SHA256)
        actual = software.get("source_sha256", {})
        controller = actual.get("batch.py")
        if controller not in (expected.get("batch.py"), PATCHED_BATCH_SHA256):
            raise ValueError("Unexpected controller fingerprint; no files changed.")
        expected["batch.py"] = controller
        if software.get("program_version") != "0.11.1" or actual != expected:
            raise ValueError("Source files differ from the reviewed v0.11.1 release; no files changed.")
    before, after = copy.deepcopy(previous), copy.deepcopy(current)
    before["source_sha256"].pop("batch.py")
    after["source_sha256"].pop("batch.py")
    if before != after:
        raise ValueError("Python, dependencies, platform or other source identity changed; use the original environment.")


def atomic_write(path, data, mode):
    """Replace one file durably without exposing a partially written document."""
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    temporary = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            os.fchmod(stream.fileno(), mode)
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        if hasattr(os, "O_DIRECTORY"):
            directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
    finally:
        temporary.unlink(missing_ok=True)


def backup(path, data, stamp):
    """Keep original bytes outside the *.py source fingerprint inventory."""
    destination = path.with_name(path.name + ".before-" + PATCH_ID + "-" + stamp + ".bak")
    if destination.exists():
        raise ValueError("Backup already exists; no file will be overwritten.")
    atomic_write(destination, data, stat.S_IMODE(path.stat().st_mode))
    return destination


def recover_batch(directory, dry_run=False):
    """Validate first, then update only the controller and its recorded identity."""
    # A helper in a subdirectory must still select the project the user entered.
    # Do not silently modify an unrelated installed copy when a source tree exists.
    project = Path.cwd()
    if (project / "polymer_clustering" / "batch.py").is_file():
        sys.path.insert(0, str(project))
    from polymer_clustering import batch
    from polymer_clustering.output import software_metadata

    directory = Path(directory).expanduser().resolve()
    manifest_path = directory / "batch.json"
    if not manifest_path.is_file():
        raise ValueError("Choose the existing interrupted batch folder containing batch.json.")
    with batch.controller_lock(directory):
        original_manifest = manifest_path.read_bytes()
        manifest = json.loads(original_manifest)
        if (manifest.get("schema_version") != 1
                or manifest.get("plan_sha256") != batch.plan_hash(manifest["jobs"])):
            raise ValueError("Batch schema or scientific plan validation failed; no files changed.")
        current = software_metadata()
        verify_software(manifest["software"], current)
        source_path = Path(batch.__file__).resolve()
        original_source = source_path.read_bytes()
        if digest(original_source) != current["source_sha256"]["batch.py"]:
            raise ValueError("Controller changed during validation; no files changed.")
        if digest(original_source) == PATCHED_BATCH_SHA256:
            updated_source = original_source
        else:
            updated_source = replacement_source(original_source)

        invocations = manifest.get("invocations")
        if not isinstance(invocations, list) or not invocations:
            raise ValueError("Expected an existing invocation record; no files changed.")
        for invocation in invocations:
            elapsed = invocation.get("wall_seconds")
            if elapsed is not None and (type(elapsed) not in (int, float)
                                       or not math.isfinite(elapsed) or elapsed < 0):
                raise ValueError("A recorded wall time is invalid; no files changed.")
        missing = sum(i.get("wall_seconds") is None for i in invocations)

        complete = 0
        for job in manifest["jobs"]:
            batch.validate_parameters(job["parameters"])
            if job["status"] == "complete":
                run = (directory / job["run_directory"]).resolve()
                if not run.is_relative_to(directory):
                    raise ValueError("Completed run points outside its batch; no files changed.")
                batch.inspect_run(run, job["parameters"], current)
                complete += 1
        print(f"Verified: {complete} completed runs; {len(manifest['jobs']) - complete} remaining.")
        print(f"Prior invocation times unavailable: {missing}; these will remain unknown.")
        print("Scientific parameters, seeds and numerical modules are unchanged.")
        if dry_run:
            print("Recovery preview passed; no source or batch record changed.")
            return
        if (original_source == updated_source
                and manifest["software"]["source_sha256"]["batch.py"] == PATCHED_BATCH_SHA256):
            print("This exact recovery fix is already applied; no changes needed.")
            return

        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + "-" + uuid.uuid4().hex[:8]
        source_backup = backup(source_path, original_source, stamp)
        manifest_backup = backup(manifest_path, original_manifest, stamp)
        updated_manifest = copy.deepcopy(manifest)
        updated_manifest["software"]["source_sha256"]["batch.py"] = PATCHED_BATCH_SHA256
        updated_manifest.setdefault("controller_patches", []).append({
            "id": PATCH_ID,
            "applied_utc": datetime.now(timezone.utc).isoformat(),
            "helper_sha256": digest(Path(__file__).read_bytes()),
            "previous_recorded_batch_sha256": manifest["software"]["source_sha256"]["batch.py"],
            "previous_installed_batch_sha256": digest(original_source),
            "patched_batch_sha256": PATCHED_BATCH_SHA256,
            "previous_last_execution_wall_seconds": manifest.get("last_execution_wall_seconds"),
            "source_backup": str(source_backup),
            "manifest_backup": str(manifest_backup),
            "scope": "Controller elapsed-time bookkeeping only; missing durations remain unknown.",
        })
        if (updated_manifest["jobs"] != manifest["jobs"]
                or updated_manifest["plan_sha256"] != manifest["plan_sha256"]):
            raise ValueError("Unexpected scientific-plan change; no files changed.")
        # Detect another writer before either replacement. A second invocation
        # can also finish a source-first patch interrupted between replacements.
        if source_path.read_bytes() != original_source or manifest_path.read_bytes() != original_manifest:
            raise ValueError("Files changed during validation; no replacements performed.")
        updated_bytes = (json.dumps(updated_manifest, indent=2, allow_nan=False) + "\n").encode()
        try:
            if updated_source != original_source:
                atomic_write(source_path, updated_source, stat.S_IMODE(source_path.stat().st_mode))
            atomic_write(manifest_path, updated_bytes, stat.S_IMODE(manifest_path.stat().st_mode))
        except BaseException:
            # For ordinary write failures restore both originals. After a hard
            # power loss the exact old/new combination is accepted on rerun.
            atomic_write(source_path, original_source, stat.S_IMODE(source_backup.stat().st_mode))
            atomic_write(manifest_path, original_manifest, stat.S_IMODE(manifest_backup.stat().st_mode))
            raise
        print("Recovery fix applied. Original controller and batch record backed up:")
        print(source_backup)
        print(manifest_backup)
        print("Use the normal --resume command in a new Python process.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", required=True, type=Path, help="Existing interrupted batch directory.")
    parser.add_argument("--dry-run", action="store_true", help="Validate without replacing source or metadata.")
    args = parser.parse_args()
    try:
        recover_batch(args.batch, args.dry_run)
    except (OSError, ValueError, KeyError, TypeError, ImportError) as error:
        print(f"Recovery stopped: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
