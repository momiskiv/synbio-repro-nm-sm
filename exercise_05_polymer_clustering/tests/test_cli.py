# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Check the runnable interface, replay, failure status, and output isolation."""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from polymer_clustering.cli import main


class CliTests(unittest.TestCase):
    def test_outputs_replay_and_no_overwrite(self) -> None:
        # Temporary directories isolate test artefacts from the user's runs.
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "input.json"
            config_path.write_text('{"seed": 1, "strain": "BB170"}', encoding="utf-8")
            output = root / "runs"
            with contextlib.redirect_stdout(io.StringIO()):
                result = main(["--config", str(config_path), "--seed", "42",
                               "--output", str(output)])
            self.assertEqual(result, 0)
            first = next(output.iterdir())
            expected = {"initial_state.npz", "bacteria.csv", "initial_state.png",
                        "initial_state.svg", "parameters.json", "run.json"}
            self.assertEqual({p.name for p in first.iterdir()}, expected)
            parameters = json.loads((first / "parameters.json").read_text())
            # The explicit flag must override the seed read from input.json.
            self.assertEqual(parameters["seed"], 42)
            metadata = json.loads((first / "run.json").read_text())
            self.assertEqual(metadata["status"], "complete")
            self.assertEqual(metadata["milestone"], "initialisation_only")
            self.assertEqual(metadata["initial_total_signals"], 0)
            timing = metadata["timing"]
            self.assertGreaterEqual(timing["total_wall_seconds"], timing["output_wall_seconds"])
            self.assertEqual(timing["simulation_wall_seconds"], 0)
            self.assertIsNone(timing["simulation_events_per_second"])
            self.assertEqual(metadata["plotting"]["requested_font_family"], "Arial")
            self.assertTrue(metadata["plotting"]["resolved_font_family"])
            # Snapshot the first run, then replay its saved settings. Numerical
            # arrays must match and none of the original files may change.
            original_files = {p.name: p.read_bytes() for p in first.iterdir()}
            with contextlib.redirect_stdout(io.StringIO()):
                result = main(["--config", str(first / "parameters.json"),
                               "--output", str(output)])
            self.assertEqual(result, 0)
            runs = list(output.iterdir())
            self.assertEqual(len(runs), 2)
            second = next(p for p in runs if p != first)
            with np.load(first / "initial_state.npz", allow_pickle=False) as a:
                with np.load(second / "initial_state.npz", allow_pickle=False) as b:
                    for key in a.files:
                        np.testing.assert_array_equal(a[key], b[key])
            for name, data in original_files.items():
                self.assertEqual((first / name).read_bytes(), data)

    def test_bad_configuration_creates_no_output(self) -> None:
        # Validation must fail before the output tree is created.
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "bad.json"
            config.write_text('{"unrecognised": 4}', encoding="utf-8")
            output = root / "outputs"
            with (
                contextlib.redirect_stderr(io.StringIO()),
                self.assertRaises(SystemExit) as error,
            ):
                main(["--config", str(config), "--output", str(output)])
            self.assertEqual(error.exception.code, 2)
            self.assertFalse(output.exists())

    def test_runtime_failure_is_recorded(self) -> None:
        # Inject a failure to check the error record without depending on a
        # particular random placement becoming jammed.
        with tempfile.TemporaryDirectory() as temporary:
            with patch(
                "polymer_clustering.cli.initialise", side_effect=ValueError("placement")
            ):
                with contextlib.redirect_stderr(io.StringIO()):
                    result = main(["--output", temporary])
            self.assertEqual(result, 1)
            directory = next(Path(temporary).iterdir())
            metadata = json.loads((directory / "run.json").read_text())
            self.assertEqual(metadata["status"], "failed")
            self.assertIn("placement", metadata["error"])
            self.assertGreaterEqual(metadata["timing"]["total_wall_seconds"], 0)
            self.assertGreaterEqual(metadata["timing"]["initialisation_wall_seconds"], 0)
            self.assertIn("finished_utc", metadata)


if __name__ == "__main__":
    unittest.main()
