# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Binding-mode controls, saved histories, replay, and incomplete-run reporting."""

import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from polymer_clustering.cli import main


class BindingCliTests(unittest.TestCase):
    def test_binding_outputs_and_config_replay(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with contextlib.redirect_stdout(io.StringIO()):
                code = main(["--config", "examples/mm32_p3_binding.json", "--duration", "0.1",
                             "--seed", "42", "--output", temporary])
            self.assertEqual(code, 0)
            first = next(Path(temporary).iterdir())
            metadata = json.loads((first / "run.json").read_text())
            self.assertEqual(metadata["milestone"], "binding_only_fixed_positions")
            self.assertEqual(metadata["binding"]["stop_time"], 0.1)
            self.assertEqual(metadata["status"], "complete")
            self.assertIn("luminescence", metadata["excluded_processes"])
            self.assertTrue({"binding_state.npz", "binding_trajectory.csv", "binding_events.csv",
                             "binding.png", "binding.svg", "polymer_binding.csv",
                             "bacteria_binding.csv"}.issubset({p.name for p in first.iterdir()}))
            settings = json.loads((first / "parameters.json").read_text())
            self.assertEqual(settings["seed"], 42)
            with (first / "binding_trajectory.csv").open() as handle:
                self.assertTrue(all(int(row["total_signals"]) == 10000
                                    for row in csv.DictReader(handle)))
            with contextlib.redirect_stdout(io.StringIO()):
                code = main(["--config", str(first / "parameters.json"), "--output", temporary])
            self.assertEqual(code, 0)
            second = next(p for p in Path(temporary).iterdir() if p != first)
            self.assertEqual((first / "binding_events.csv").read_bytes(),
                             (second / "binding_events.csv").read_bytes())
            with np.load(first / "binding_state.npz", allow_pickle=False) as a:
                with np.load(second / "binding_state.npz", allow_pickle=False) as b:
                    for key in a.files:
                        np.testing.assert_array_equal(a[key], b[key])

    def test_incomplete_run_has_nonzero_exit_and_no_future_observations(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with contextlib.redirect_stdout(io.StringIO()):
                code = main(["--mode", "binding", "--duration", "100", "--max-events", "1",
                             "--output", temporary])
            self.assertEqual(code, 3)
            directory = next(Path(temporary).iterdir())
            metadata = json.loads((directory / "run.json").read_text())
            self.assertEqual(metadata["status"], "incomplete")
            self.assertEqual(metadata["binding"]["status"], "event_limit")
            with (directory / "binding_trajectory.csv").open() as handle:
                times = [float(row["time"]) for row in csv.DictReader(handle)]
            self.assertEqual(times[-1], metadata["binding"]["stop_time"])
            self.assertLess(times[-1], 100)

    def test_invalid_binding_controls_create_no_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            for options in (["--mode", "binding", "--polymers", "10"],
                            ["--mode", "binding", "--k-bs-on", "nan"],
                            ["--duration", "10"]):
                output = Path(temporary) / "must-not-exist"
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    main([*options, "--output", str(output)])
                self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
