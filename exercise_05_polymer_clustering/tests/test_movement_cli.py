# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Movement-mode output, settings replay, compatibility, and status checks."""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from PIL import Image

from polymer_clustering.cli import main


class MovementCliTests(unittest.TestCase):
    def test_output_animation_and_numerical_replay(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            arguments = ["--config", "examples/mm32_p3_movement.json", "--duration", "0.05",
                         "--spatial-frames", "3", "--output", temporary]
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(arguments), 0)
            first = next(Path(temporary).iterdir())
            metadata = json.loads((first / "run.json").read_text())
            self.assertEqual(metadata["milestone"], "binding_and_movement")
            self.assertNotIn("movement", metadata["excluded_processes"])
            self.assertGreater(metadata["movement_event_counts"]["signal_hop"], 0)
            self.assertGreater(metadata["binding"]["events"], metadata["binding"]["recorded_events"])
            expected = {"movement.png", "movement.svg", "movement.gif", "movement_positions.csv",
                        "movement_frames.npz", "final_bacteria.csv"}
            self.assertTrue(expected.issubset({p.name for p in first.iterdir()}))
            with Image.open(first / "movement.gif") as gif:
                self.assertEqual(gif.n_frames, 3)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--config", str(first / "parameters.json"),
                                       "--output", temporary]), 0)
            second = next(p for p in Path(temporary).iterdir() if p != first)
            with np.load(first / "binding_state.npz", allow_pickle=False) as a:
                with np.load(second / "binding_state.npz", allow_pickle=False) as b:
                    for key in a.files:
                        np.testing.assert_array_equal(a[key], b[key])

    def test_legacy_file_missing_rng_scheme_and_invalid_mode_controls(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "v02.json"
            config.write_text('{"mode":"initial","seed":2026}')
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--config", str(config), "--output", str(root / "runs")]), 0)
            run = next((root / "runs").iterdir())
            self.assertEqual(json.loads((run / "parameters.json").read_text())["rng_scheme"], "legacy_v02")
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(["--mode", "binding", "--d-s", "1", "--output", str(root / "bad")])
            self.assertFalse((root / "bad").exists())

    def test_incomplete_movement_run_is_saved(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--mode", "movement", "--max-events", "1",
                                       "--spatial-frames", "2", "--output", temporary]), 3)
            directory = next(Path(temporary).iterdir())
            metadata = json.loads((directory / "run.json").read_text())
            self.assertEqual(metadata["status"], "incomplete")
            self.assertEqual(metadata["binding"]["events"], 1)
            with np.load(directory / "movement_frames.npz", allow_pickle=False) as frames:
                self.assertEqual(frames["times"][-1], metadata["binding"]["stop_time"])


if __name__ == "__main__":
    unittest.main()
