# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Regression checks using balanced fixtures; no assignment trajectories run."""

import contextlib
import copy
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
PROJECT = Path(os.environ.get("V0111_PROJECT", Path.cwd())).resolve()
sys.path.insert(0, str(PROJECT))
sys.path.insert(0, str(HERE))
import apply_resume_timing_fix as recovery
from polymer_clustering import batch
from polymer_clustering.batch_output import save_ensemble
from polymer_clustering.output import software_metadata, write_json
from tests.test_batch import small_spec, fixture_run


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.old_source = (PROJECT / "polymer_clustering/batch.py").read_bytes()
        self.patched_source = recovery.replacement_source(self.old_source)
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "module/batch.py"
        self.source.parent.mkdir()
        self.source.write_bytes(self.old_source)
        self.directory = self.root / "batch"
        self.directory.mkdir()
        self.software = software_metadata()
        self.spec = small_spec()
        self.jobs = batch.make_plan(self.spec)
        # Two completed jobs and two unfinished jobs exercise retained results.
        for job in self.jobs[:2]:
            path = fixture_run(self.directory, job, self.software)
            job.update(status="complete", run_directory=str(path.relative_to(self.directory)))
        self.jobs[2].update(status="running", attempts=[{"directory": "old_attempt"}])
        self.manifest = dict(schema_version=1, status="running", jobs=self.jobs,
                             plan_sha256=batch.plan_hash(self.jobs), software=self.software,
                             invocations=[dict(started_utc="2026-10-22T12:27:01+00:00", workers=8)])
        self.manifest_path = self.directory / "batch.json"
        write_json(self.manifest_path, self.manifest)

    def metadata(self):
        data = copy.deepcopy(self.software)
        data["source_sha256"]["batch.py"] = recovery.digest(self.source.read_bytes())
        return data

    def apply(self, preview=False):
        with patch.object(batch, "__file__", str(self.source)), \
             patch("polymer_clustering.output.software_metadata", side_effect=self.metadata), \
             contextlib.redirect_stdout(io.StringIO()):
            recovery.recover_batch(self.directory, preview)

    def module(self):
        path = self.root / "patched.py"
        path.write_bytes(self.patched_source)
        spec = importlib.util.spec_from_file_location("polymer_clustering._recovery_test", path)
        result = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(result)
        return result

    def test_preview_is_read_only(self):
        old_json = self.manifest_path.read_bytes()
        self.apply(preview=True)
        self.assertEqual(self.manifest_path.read_bytes(), old_json)
        self.assertEqual(self.source.read_bytes(), self.old_source)
        self.assertFalse(list(self.root.rglob("*.bak")))

    def test_install_preserves_jobs_results_and_unknown_time_and_is_idempotent(self):
        old_json = self.manifest_path.read_bytes()
        original_runs = {p: p.read_bytes() for p in self.directory.rglob("*.csv")}
        self.apply()
        data = json.loads(self.manifest_path.read_text())
        self.assertEqual(data["jobs"], self.manifest["jobs"])
        self.assertEqual(data["plan_sha256"], self.manifest["plan_sha256"])
        self.assertEqual(data["invocations"], self.manifest["invocations"])
        self.assertEqual(self.source.read_bytes(), self.patched_source)
        self.assertEqual(data["software"]["source_sha256"]["batch.py"], recovery.PATCHED_BATCH_SHA256)
        for name in batch.NUMERICAL_MODULES:
            self.assertEqual(data["software"]["source_sha256"][name], self.software["source_sha256"][name])
        for path, original in original_runs.items():
            self.assertEqual(path.read_bytes(), original)
        receipt = data["controller_patches"][-1]
        self.assertEqual(Path(receipt["source_backup"]).read_bytes(), self.old_source)
        self.assertEqual(Path(receipt["manifest_backup"]).read_bytes(), old_json)
        fixed_json = self.manifest_path.read_bytes()
        self.apply()
        self.assertEqual(self.manifest_path.read_bytes(), fixed_json)
        self.assertEqual(len(list(self.root.rglob("*.bak"))), 2)

    def test_unsupported_source_and_changed_environment_rejected(self):
        changed = self.old_source + b"\n# unexpected edit\n"
        self.source.write_bytes(changed)
        with self.assertRaises(ValueError):
            self.apply()
        self.assertEqual(self.source.read_bytes(), changed)
        self.source.write_bytes(self.old_source)
        data = copy.deepcopy(self.manifest)
        data["software"]["numpy"] = "different"
        write_json(self.manifest_path, data)
        with self.assertRaises(ValueError):
            self.apply()
        self.assertEqual(self.source.read_bytes(), self.old_source)

    def test_changed_plan_and_corrupt_completed_output_rejected(self):
        data = copy.deepcopy(self.manifest)
        data["jobs"][0]["parameters"]["seed"] += 1
        write_json(self.manifest_path, data)
        with self.assertRaises(ValueError):
            self.apply()
        write_json(self.manifest_path, self.manifest)
        run = self.directory / self.jobs[0]["run_directory"] / "run.json"
        run.write_text("{}")
        with self.assertRaises(ValueError):
            self.apply()
        self.assertEqual(self.source.read_bytes(), self.old_source)

    def test_active_controller_rejected(self):
        with batch.controller_lock(self.directory), self.assertRaises(ValueError):
            self.apply()
        self.assertEqual(self.source.read_bytes(), self.old_source)

    def test_source_first_power_loss_can_be_repaired(self):
        self.source.write_bytes(self.patched_source)
        self.apply()
        data = json.loads(self.manifest_path.read_text())
        self.assertEqual(data["jobs"], self.manifest["jobs"])
        self.assertEqual(data["software"]["source_sha256"]["batch.py"], recovery.PATCHED_BATCH_SHA256)

    def test_write_failure_rolls_back_both_files(self):
        original_json = self.manifest_path.read_bytes()
        actual_write = recovery.atomic_write
        failed = False
        def fail_once(path, data, mode):
            nonlocal failed
            if path == self.manifest_path and not failed:
                failed = True
                raise OSError("simulated write failure")
            return actual_write(path, data, mode)
        with patch.object(recovery, "atomic_write", side_effect=fail_once), self.assertRaises(OSError):
            self.apply()
        self.assertEqual(self.source.read_bytes(), self.old_source)
        self.assertEqual(self.manifest_path.read_bytes(), original_json)

    def test_unknown_duration_stays_unknown(self):
        module = self.module()
        data = {"invocations": [{"started_utc": "earlier"}, {"wall_seconds": 12.5}]}
        module.record_controller_wall_time(data)
        self.assertNotIn("wall_seconds", data["invocations"][0])
        self.assertEqual(data["invocations"][0]["timing_status"], "unavailable_after_interruption")
        self.assertIsNone(data["controller_wall_seconds_all_invocations"])
        self.assertEqual(data["controller_wall_seconds_recorded_invocations"], 12.5)
        self.assertEqual(data["controller_wall_seconds_missing_invocations"], 1)

    def test_known_durations_and_invalid_durations(self):
        module = self.module()
        data = {"invocations": [{"wall_seconds": 1.25}, {"wall_seconds": 2.75}]}
        module.record_controller_wall_time(data)
        self.assertEqual(data["controller_wall_seconds_all_invocations"], 4)
        self.assertEqual(data["controller_wall_seconds_missing_invocations"], 0)
        for invalid in (True, -1, float("nan"), float("inf"), "12"):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                module.record_controller_wall_time({"invocations": [{"wall_seconds": invalid}]})

    def test_resumed_batch_finalizes_ensemble_with_missing_old_duration(self):
        module = self.module()
        # Make the remaining balanced fixture runs complete without SSA execution.
        for job in self.jobs[2:]:
            path = fixture_run(self.directory, job, self.software)
            job.update(status="complete", run_directory=str(path.relative_to(self.directory)))
        self.manifest["software"] = copy.deepcopy(self.software)
        self.manifest["software"]["source_sha256"]["batch.py"] = recovery.PATCHED_BATCH_SHA256
        with patch.object(module, "execute_jobs", return_value=0), \
             contextlib.redirect_stdout(io.StringIO()):
            result = module.finish_batch(self.directory, self.manifest, [], 4, time.perf_counter(), "off")
        self.assertEqual(result, 0)
        saved = json.loads(self.manifest_path.read_text())
        self.assertEqual(saved["status"], "complete")
        self.assertIsNone(saved["controller_wall_seconds_all_invocations"])
        self.assertGreaterEqual(saved["invocations"][-1]["wall_seconds"], 0)
        self.assertEqual(saved["controller_wall_seconds_missing_invocations"], 1)
        self.assertTrue((self.directory / "ensemble_trajectory.csv").is_file())
        self.assertTrue((self.directory / "ensemble.png").is_file())


if __name__ == "__main__":
    unittest.main()
