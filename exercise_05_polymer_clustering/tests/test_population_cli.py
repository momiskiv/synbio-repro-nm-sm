# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Full-mode files, changing-population archives/GIFs, replay, and invalid controls."""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from PIL import Image

from polymer_clustering.cli import main


class PopulationCliTests(unittest.TestCase):
    def test_full_mode_birth_frames_outputs_and_replay(self):
        with tempfile.TemporaryDirectory() as temporary:
            options = ["--mode", "full", "--grid-size", "8", "--bacteria", "1",
                       "--strain", "BB170", "--duration", "4", "--sample-interval", "1",
                       "--doubling-time", "2", "--division-width", "0", "--population-cap", "4",
                       "--production-rate", "2", "--spatial-frames", "3", "--output", temporary]
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(options), 0)
            run = next(Path(temporary).iterdir())
            metadata = json.loads((run / "run.json").read_text())
            self.assertEqual(metadata["status"], "complete")
            self.assertEqual(metadata["milestone"], "population_coupled_model")
            self.assertEqual(metadata["excluded_processes"], [])
            self.assertEqual(metadata["population_model"]["division_support"], [2, 2])
            self.assertEqual(metadata["population"]["bacteria"], 4)
            self.assertEqual(metadata["population"]["actual_total_signals"],
                             metadata["population"]["produced_signals"])
            expected = {"population_trajectory.csv", "lineage.csv", "population.png", "population.svg",
                        "response.gif", "movement.gif"}
            self.assertTrue(expected <= {p.name for p in run.iterdir()})
            with np.load(run / "movement_frames.npz", allow_pickle=False) as frames:
                np.testing.assert_array_equal(frames["bacteria_counts"], [1, 2, 4])
                self.assertEqual(frames["bacteria_positions"].shape, (3, 4, 3))
                self.assertTrue(np.all(frames["bacteria_positions"][0, 1:] == -1))
            with np.load(run / "response_frames.npz", allow_pickle=False) as frames:
                self.assertTrue(np.all(np.isnan(frames["expression"][0, 1:])))
                self.assertTrue(np.all(np.isfinite(frames["luminescence"][-1])))
                np.testing.assert_array_equal(frames["present"].sum(axis=1), [1, 2, 4])
            for name in ("response.gif", "movement.gif"):
                with Image.open(run / name) as gif:
                    self.assertEqual(gif.n_frames, 3)
                    dimensions = gif.size
                    for index in range(3):
                        gif.seek(index)
                        self.assertEqual(gif.size, dimensions)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--config", str(run / "parameters.json"),
                                       "--output", temporary]), 0)
            replay = next(p for p in Path(temporary).iterdir() if p != run)
            for filename in ("binding_events.csv", "population_trajectory.csv", "lineage.csv"):
                self.assertEqual((run / filename).read_bytes(), (replay / filename).read_bytes())

    def test_invalid_population_settings_leave_no_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "bad"
            for flags in (["--mode", "response", "--no-growth"],
                          ["--mode", "full", "--bacteria", "25", "--population-cap", "10"],
                          ["--mode", "full", "--division-width", "30002"]):
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    main([*flags, "--output", str(output)])
                self.assertFalse(output.exists())

    def test_full_mode_event_cap_is_incomplete(self):
        with tempfile.TemporaryDirectory() as temporary:
            with contextlib.redirect_stdout(io.StringIO()):
                code = main(["--mode", "full", "--grid-size", "8", "--bacteria", "1",
                             "--strain", "BB170", "--duration", "10", "--production-rate", "0",
                             "--d-b", "0", "--doubling-time", "2", "--division-width", "0",
                             "--max-events", "1", "--spatial-frames", "2", "--output", temporary])
            self.assertEqual(code, 3)
            metadata = json.loads((next(Path(temporary).iterdir()) / "run.json").read_text())
            self.assertEqual(metadata["status"], "incomplete")
            self.assertEqual(metadata["population"]["time"], 2)


if __name__ == "__main__":
    unittest.main()
