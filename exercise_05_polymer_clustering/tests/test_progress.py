# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Check progress output, failures, and exact seeded numerical equivalence."""

import contextlib
import copy
from dataclasses import fields, is_dataclass
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from polymer_clustering.batch import (
    PARAMETERISED_BINDING_SHA256, PRE_PROGRESS_BINDING_SHA256,
    PRE_PROGRESS_MOVEMENT_SHA256, OPTIMISED_MOVEMENT_SHA256, available_cpus,
    main as batch_main, numerical_sources_match,
)
from polymer_clustering.binding import BindingConfig, simulate_binding
from polymer_clustering.cli import main
from polymer_clustering.model import InitialConfig, initialise
from polymer_clustering.movement import MovementConfig, simulate_movement
from polymer_clustering.output import software_metadata, write_json
from polymer_clustering.population import PopulationConfig
from polymer_clustering.progress import BatchProgress, ConsoleProgress, RunProgress
from polymer_clustering.response import ResponseConfig
from test_batch import small_spec


class TTY(io.StringIO):
    def isatty(self):
        return True


class ProgressChecks(unittest.TestCase):
    def assert_same(self, left, right):
        """Compare every result field, including event history and per-cell arrays."""
        if isinstance(left, np.ndarray):
            np.testing.assert_array_equal(left, right)
        elif is_dataclass(left):
            for field in fields(left):
                self.assert_same(getattr(left, field.name), getattr(right, field.name))
        elif isinstance(left, dict):
            self.assertEqual(left.keys(), right.keys())
            for key in left:
                self.assert_same(left[key], right[key])
        elif isinstance(left, (list, tuple)):
            self.assertEqual(len(left), len(right))
            for a, b in zip(left, right):
                self.assert_same(a, b)
        else:
            self.assertEqual(left, right)

    def test_tty_redraw_is_throttled_bounded_and_ends_with_newline(self):
        stream, now = TTY(), [0.]
        progress = ConsoleProgress("auto", stream, lambda: now[0])
        with patch("polymer_clustering.progress.shutil.get_terminal_size",
                   return_value=os.terminal_size((40, 24))):
            progress.update("a" * 100)
            progress.update("not yet")
            now[0] = 1
            progress.update("short")
            progress.close()
        value = stream.getvalue()
        self.assertNotIn("not yet", value)
        self.assertNotIn("\x1b", value)
        self.assertTrue(value.endswith("\n"))
        self.assertTrue(all(len(line) <= 39 for line in value.rstrip("\n").split("\r")))

    def test_non_tty_modes_and_plain_log_throttling(self):
        for mode in ("auto", "off"):
            stream = io.StringIO()
            progress = RunProgress(10, mode, stream=stream)
            self.assertIsNone(progress.callback)
            progress.phase("simulating")
            progress.finish("complete")
            self.assertEqual(stream.getvalue(), "")
        stream, now = io.StringIO(), [0.]
        progress = ConsoleProgress("on", stream, lambda: now[0])
        progress.update("first")
        now[0] = 5
        progress.update("throttled")
        now[0] = 10
        progress.update("later")
        progress.close()
        self.assertEqual(stream.getvalue(), "first\nlater\n")

    def test_worker_status_is_throttled_atomic_and_marks_actual_incomplete_time(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, now = Path(tmp)/"progress.json", [0.]
            progress = RunProgress(10, "off", path, clock=lambda: now[0])
            self.assertIsNotNone(progress.callback)
            progress.phase("simulating")
            progress.observe(4, 1000)
            self.assertEqual(json.loads(path.read_text())["model_time"], 0)
            now[0] = 1
            progress.observe(5, 2000)
            self.assertEqual(json.loads(path.read_text())["model_time"], 5)
            progress.finish("incomplete")
            record = json.loads(path.read_text())
            self.assertEqual(record["phase"], "incomplete")
            self.assertEqual(record["model_time"], 5)
            self.assertEqual(record["events"], 2000)
            self.assertFalse(path.with_suffix(".json.tmp").exists())

    def test_progress_io_failures_do_not_raise_into_the_solver(self):
        class ClosedStream(io.StringIO):
            def write(self, value):
                raise OSError("closed terminal")
        progress = RunProgress(10, "on", stream=ClosedStream())
        progress.phase("simulating")
        self.assertFalse(progress.console.enabled)
        with tempfile.TemporaryDirectory() as tmp:
            stream = io.StringIO()
            progress = RunProgress(10, "off", Path(tmp)/"missing"/"progress.json", stream=stream)
            progress.observe(1, 1000)
            progress.finish("complete")
            self.assertIsNone(progress.status_file)
            self.assertIn("Progress file unavailable", stream.getvalue())

    def test_batch_counts_verified_runs_and_accepts_missing_or_invalid_worker_status(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()) as out:
            root = Path(tmp)
            jobs = [dict(id="reuse/seed_1", status="complete", attempts=[]),
                    dict(id="new/seed_1", status="running", attempts=[dict(directory="attempt")])]
            display = BatchProgress(root, jobs, "on")
            display.update(force=True)
            attempt = root/"attempt"
            attempt.mkdir()
            (attempt/"progress.json").write_text("not valid JSON")
            display.update(force=True)
            write_json(attempt/"progress.json", dict(phase="saving outputs", model_time=10, duration=10))
            display.update(force=True)
            jobs[1]["status"] = "incomplete"
            display.close()
            value = out.getvalue()
            self.assertIn("1/2 complete", value)
            self.assertIn("saving outputs", value)
            self.assertIn("1 incomplete", value)
            self.assertNotIn("2/2 complete", value)
            self.assertNotIn("\r", value)

    def test_first_worker_percentage_is_visible_in_an_eighty_column_terminal(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(TTY()) as out:
            root = Path(tmp)
            (root/"attempt").mkdir()
            write_json(root/"attempt"/"progress.json",
                       dict(phase="simulating", model_time=2500, duration=30000))
            jobs = [dict(id="p3_3000/seed_2026", status="running",
                         attempts=[dict(directory="attempt")])]
            display = BatchProgress(root, jobs, "auto")
            with patch("polymer_clustering.progress.shutil.get_terminal_size",
                       return_value=os.terminal_size((80, 24))):
                display.update(force=True)
                display.close()
            self.assertIn("simulation 8.3%", out.getvalue())
            self.assertTrue(out.getvalue().endswith("\n"))

    def test_all_solver_fields_match_with_read_only_progress_for_dense_and_capped_runs(self):
        for cutoff, cap in ((10, 100000), (None, 100000), (10, 1)):
            initial = initialise(InitialConfig(grid_size=6, bacteria=2, polymers=360,
                                               seed=17, rng_scheme="independent_v1"))
            binding = BindingConfig(duration=2, sample_interval=.5, max_events=cap,
                                    k_pb_on=.005, k_pb_off=.001, k_ps_on=.01, k_ps_off=.001,
                                    ps_crowding_threshold=cutoff)
            for moving in (False, True):
                with self.subTest(cutoff=cutoff, cap=cap, moving=moving):
                    calls = []
                    def run(callback):
                        if moving:
                            return simulate_movement(initial, binding,
                                MovementConfig(spatial_frames=3, event_log="all"), 17,
                                response=ResponseConfig(),
                                population=PopulationConfig(doubling_time=1, division_width=0,
                                                            population_cap=4, production_rate=10),
                                progress=callback)
                        return simulate_binding(initial, binding, 17, progress=callback)
                    plain = run(None)
                    observed = run(lambda t, n: calls.append((t, n)))
                    self.assert_same(plain, observed)
                    count = observed.total_events if moving else len(observed.events)
                    self.assertEqual(calls[0], (0., 0))
                    self.assertEqual(calls[-1], (observed.time, count))
                    self.assertTrue(all(n % 1000 == 0 for _, n in calls[1:-1]))
                    self.assertEqual(observed.status, "event_limit" if cap == 1 else "complete")

    def test_exact_progress_fingerprints_reuse_both_v09_cutoffs_and_reject_edits(self):
        current = software_metadata()
        self.assertEqual(current["source_sha256"]["binding.py"], PARAMETERISED_BINDING_SHA256)
        self.assertEqual(current["source_sha256"]["movement.py"], OPTIMISED_MOVEMENT_SHA256)
        old = copy.deepcopy(current)
        old["source_sha256"].update({"binding.py": PRE_PROGRESS_BINDING_SHA256,
                                     "movement.py": PRE_PROGRESS_MOVEMENT_SHA256})
        for cutoff in (10, None, 20):
            self.assertTrue(numerical_sources_match(old, current, {"ps_crowding_threshold": cutoff}))
        for name in ("binding.py", "movement.py", "model.py"):
            changed = copy.deepcopy(old)
            changed["source_sha256"][name] = "unreviewed edit"
            self.assertFalse(numerical_sources_match(changed, current, {}))
            changed = copy.deepcopy(current)
            changed["source_sha256"][name] = "unreviewed target"
            self.assertFalse(numerical_sources_match(old, changed, {}))

    def test_cli_event_limit_failure_interrupt_and_non_model_progress_setting(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = root/"config.json"
            settings = dict(mode="binding", grid_size=4, bacteria=1, duration=100,
                            sample_interval=1, max_events=1, font_family="DejaVu Sans")
            write_json(config, settings)
            status = root/"progress.json"
            for name, error, expected in (("limited", None, 3), ("failed", ValueError("check"), 1),
                                           ("interrupted", KeyboardInterrupt(), 130)):
                with contextlib.ExitStack() as stack:
                    stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                    stack.enter_context(contextlib.redirect_stderr(io.StringIO()))
                    if error is not None:
                        stack.enter_context(patch("polymer_clustering.cli.simulate_binding",
                                                  side_effect=error))
                    code = main(["--config", str(config), "--output", str(root/name),
                                 "--progress", "on", "--progress-file", str(status)])
                self.assertEqual(code, expected)
                meta = json.loads(next((root/name).glob("*/run.json")).read_text())
                record = json.loads(status.read_text())
                self.assertEqual(record["phase"], meta["status"])
                self.assertNotEqual(record["phase"], "complete")
                self.assertNotIn("progress", meta["parameters"])
                if name == "limited":
                    self.assertLess(record["model_time"], record["duration"])
            write_json(config, {**settings, "progress": "on"})
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(["--config", str(config), "--output", str(root/"invalid")])
            self.assertFalse((root/"invalid").exists())

    def test_actual_workers_on_and_off_match_and_resume_does_not_relaunch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            spec = small_spec()
            spec["seeds"] = [21]
            config = root/"batch.json"
            write_json(config, spec)
            batches = []
            for mode in ("off", "on"):
                with contextlib.redirect_stdout(io.StringIO()) as output:
                    self.assertEqual(batch_main(["--config", str(config), "--progress", mode,
                        "--workers", str(min(2, available_cpus())), "--output", str(root/mode)]), 0)
                directory = next((root/mode).iterdir())
                manifest = json.loads((directory/"batch.json").read_text())
                batches.append((directory, manifest))
                for job in manifest["jobs"]:
                    status = directory/job["attempts"][0]["directory"]/"progress.json"
                    if mode == "on":
                        self.assertEqual(json.loads(status.read_text())["phase"], "complete")
                    else:
                        self.assertFalse(status.exists())
                if mode == "on":
                    self.assertIn("2/2 complete", output.getvalue())
                    self.assertNotIn("\r", output.getvalue())
            for left, right in zip(batches[0][1]["jobs"], batches[1][1]["jobs"]):
                a = batches[0][0]/left["run_directory"]
                b = batches[1][0]/right["run_directory"]
                for path in a.glob("*.csv"):
                    self.assertEqual(path.read_bytes(), (b/path.name).read_bytes(), path.name)
                for path in a.glob("*.npz"):
                    with np.load(path, allow_pickle=False) as aa, np.load(b/path.name, allow_pickle=False) as bb:
                        self.assertEqual(aa.files, bb.files)
                        for key in aa.files:
                            np.testing.assert_array_equal(aa[key], bb[key])
            directory, manifest = batches[1]
            with patch("polymer_clustering.batch.subprocess.Popen") as launch, \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(batch_main(["--resume", str(directory), "--progress", "off"]), 0)
                launch.assert_not_called()
            after = json.loads((directory/"batch.json").read_text())
            self.assertEqual([j["attempts"] for j in manifest["jobs"]],
                             [j["attempts"] for j in after["jobs"]])


if __name__ == "__main__":
    unittest.main()
