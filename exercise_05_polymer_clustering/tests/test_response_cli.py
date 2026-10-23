# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Response controls, outputs, replay, and incomplete-run reporting."""

import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from PIL import Image

from polymer_clustering.cli import main


class ResponseCliTests(unittest.TestCase):
    def test_response_outputs_metadata_and_replay(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with contextlib.redirect_stdout(io.StringIO()):
                code = main(["--config", "examples/mm32_p3_response.json", "--duration", "0.1",
                             "--spatial-frames", "3", "--light-map-max", "5000",
                             "--output", temporary])
            self.assertEqual(code, 0)
            first = next(Path(temporary).iterdir())
            metadata = json.loads((first / "run.json").read_text())
            self.assertEqual(metadata["milestone"], "binding_movement_and_response")
            self.assertEqual(metadata["excluded_processes"], ["growth", "signal_production"])
            self.assertEqual(metadata["response_model"]["theoretical_light_upper_limit_per_cell"], 5000)
            self.assertEqual(metadata["response"]["time"], 0.1)
            self.assertEqual(metadata["response_model"]["light_map_max"], 5000)
            expected = {"response_state.npz", "response_frames.npz", "cell_response.csv",
                        "response_trajectory.csv", "response.png", "response.svg",
                        "response_spatial.png", "response_spatial.svg", "response.gif"}
            self.assertTrue(expected.issubset({p.name for p in first.iterdir()}))
            with Image.open(first / "response.gif") as gif:
                self.assertEqual(gif.n_frames, 3)
            with (first / "cell_response.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 50)  # t=0 and the final, off-grid observation.
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--config", str(first / "parameters.json"),
                                       "--output", temporary]), 0)
            second = next(p for p in Path(temporary).iterdir() if p != first)
            with np.load(first / "response_state.npz", allow_pickle=False) as a:
                with np.load(second / "response_state.npz", allow_pickle=False) as b:
                    for key in a.files:
                        np.testing.assert_array_equal(a[key], b[key])

    def test_invalid_response_controls_create_no_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "bad"
            for options in (["--mode", "response", "--k-deg", "-0.01"],
                            ["--mode", "movement", "--k-syn", "1"]):
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    main([*options, "--output", str(output)])
                self.assertFalse(output.exists())

    def test_incomplete_response_records_actual_stop(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with contextlib.redirect_stdout(io.StringIO()):
                code = main(["--mode", "response", "--max-events", "1",
                             "--spatial-frames", "2", "--output", temporary])
            self.assertEqual(code, 3)
            run = next(Path(temporary).iterdir())
            metadata = json.loads((run / "run.json").read_text())
            self.assertEqual(metadata["status"], "incomplete")
            self.assertGreater(metadata["timing"]["simulation_wall_seconds"], 0)
            self.assertGreater(metadata["timing"]["simulation_events_per_second"], 0)
            self.assertGreaterEqual(metadata["timing"]["total_wall_seconds"],
                                    metadata["timing"]["simulation_wall_seconds"] +
                                    metadata["timing"]["output_wall_seconds"])
            self.assertEqual(metadata["response"]["time"], metadata["binding"]["stop_time"])


if __name__ == "__main__":
    unittest.main()
