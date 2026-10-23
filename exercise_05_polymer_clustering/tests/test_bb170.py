# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Check BB170 source conditions, lifetime boundaries and exact old-run reuse."""

import copy
import contextlib
import io
import json
from pathlib import Path
import unittest
import tempfile

import numpy as np

from polymer_clustering.batch import (
    make_plan, numerical_sources_match, PRE_BOUNDARY_POPULATION_SHA256,
    POSITIVE_LIFETIME_POPULATION_SHA256,
)
from polymer_clustering.binding import BindingConfig
from polymer_clustering.cli import main
from polymer_clustering.batch_output import save_ensemble
from polymer_clustering.model import InitialConfig, initialise
from polymer_clustering.movement import MovementConfig, simulate_movement
from polymer_clustering.output import software_metadata
from polymer_clustering.population import PopulationConfig
from polymer_clustering.response import ResponseConfig
from test_batch import fixture_run, small_spec


ROOT = Path(__file__).resolve().parents[1]


class BB170Checks(unittest.TestCase):
    def test_existing_positive_support_lifetimes_preserve_values_and_rng_state(self):
        # Freeze the predecessor's direct NumPy sampling as an independent oracle.
        # This covers defaults, a different valid width, and deterministic times.
        for mean, width in [(15000., 20000.), (10000., 4000.), (2., 0.)]:
            for seed in [17, 2026, 2035]:
                with self.subTest(mean=mean, width=width, seed=seed):
                    actual_rng = np.random.default_rng(seed)
                    old_rng = np.random.default_rng(seed)
                    config = PopulationConfig(doubling_time=mean, division_width=width)
                    actual = [config.lifetime(actual_rng) for _ in range(100)]
                    expected = ([float(old_rng.triangular(mean-width/2, mean, mean+width/2))
                                 for _ in range(100)] if width else [mean] * 100)
                    np.testing.assert_array_equal(actual, expected)
                    self.assertEqual(actual_rng.bit_generator.state, old_rng.bit_generator.state)

    def test_zero_support_boundary_matches_continuous_triangular_distribution(self):
        config = PopulationConfig(doubling_time=10000, division_width=20000)
        rng = np.random.default_rng(32)
        times = np.array([config.lifetime(rng) for _ in range(20000)])
        self.assertGreater(times.min(), 0)
        self.assertLessEqual(times.max(), 20000)
        self.assertAlmostEqual(times.mean(), 10000, delta=100)
        self.assertAlmostEqual(times.var(), 20000**2 / 24, delta=600000)
        for parameters in [dict(doubling_time=9999, division_width=20000),
                           dict(doubling_time=0, division_width=0)]:
            with self.assertRaises(ValueError):
                PopulationConfig(**parameters)

    def test_finite_precision_zero_is_redrawn_before_scheduling(self):
        class Draw:
            def __init__(self):
                self.values = iter([0., 4000.])
                self.bounds = []

            def triangular(self, left, mode, right):
                self.bounds.append((left, mode, right))
                return next(self.values)

        draw = Draw()
        self.assertEqual(PopulationConfig(doubling_time=10000,
                                         division_width=20000).lifetime(draw), 4000)
        self.assertEqual(draw.bounds, [(0., 10000, 20000.), (0., 10000, 20000.)])

    def test_boundary_growth_runs_with_production_and_conservation(self):
        # Accelerated timing is solely a scheduling diagnostic, not scientific data.
        initial = initialise(InitialConfig(grid_size=8, bacteria=1, polymers=3,
                                           strain="BB170", seed=2026))
        result = simulate_movement(
            initial, BindingConfig(duration=6, sample_interval=.5,
                                   k_pb_on=.005, k_pb_off=.001,
                                   k_ps_on=.01, k_ps_off=.001),
            MovementConfig(d_s=.1, d_b=.05, validate_every=1, spatial_frames=5),
            7, ResponseConfig(), "BB170",
            PopulationConfig(doubling_time=1, division_width=2,
                             population_cap=4, production_rate=2),
        )
        self.assertEqual(result.status, "complete")
        self.assertEqual(len(result.state.bacteria), 4)
        self.assertGreater(result.produced_signals, 0)
        self.assertTrue(all(row["birth_time"] > 0 for row in result.lineage
                            if row["parent"] >= 0))
        self.assertTrue(all(row["actual_total_signals"] == row["expected_total_signals"]
                            for row in result.population_trajectory))
        result.state.validate(result.produced_signals, 3)

    def test_old_population_fingerprint_bridge_is_exact_and_parameter_restricted(self):
        current = software_metadata()
        self.assertEqual(current["source_sha256"]["population.py"],
                         POSITIVE_LIFETIME_POPULATION_SHA256)
        old = copy.deepcopy(current)
        old["source_sha256"]["population.py"] = PRE_BOUNDARY_POPULATION_SHA256
        for parameters in [{}, dict(doubling_time=15000, division_width=20000),
                           dict(doubling_time=2, division_width=0)]:
            self.assertTrue(numerical_sources_match(old, current, parameters))
        for mean in [10000, 9999]:
            self.assertFalse(numerical_sources_match(
                old, current, dict(doubling_time=mean, division_width=20000)))
        edited = copy.deepcopy(old)
        edited["source_sha256"]["population.py"] = "unreviewed predecessor"
        self.assertFalse(numerical_sources_match(edited, current, {}))
        edited = copy.deepcopy(current)
        edited["source_sha256"]["population.py"] = "unreviewed target"
        self.assertFalse(numerical_sources_match(old, edited, {}))
        old["numpy"] = "another NumPy version"
        self.assertFalse(numerical_sources_match(old, current, {}))

    def test_collection_plans_cover_source_loadings_and_ten_seeds(self):
        for regime, bacteria, mean in [("standard", 25, 15000),
                                       ("low_density", 2, 15000),
                                       ("faster_growth", 25, 10000)]:
            spec = json.loads((ROOT/f"examples/bb170_{regime}_full.json").read_text())
            jobs = make_plan(spec)
            self.assertEqual(len(jobs), 80)
            self.assertEqual(len({j["id"] for j in jobs}), 80)
            self.assertEqual(spec["seeds"], list(range(2026, 2036)))
            self.assertEqual({j["parameters"]["polymers"] for j in jobs}, set(range(0, 7001, 1000)))
            for job in jobs:
                p = job["parameters"]
                self.assertEqual((p["strain"], p["bacteria"], p["doubling_time"],
                                  p["division_width"]), ("BB170", bacteria, mean, 20000))
                self.assertEqual((p["duration"], p["sample_interval"], p["spatial_frames"]),
                                 (30000, 100, 61))
                self.assertEqual((p["k_pb_on"], p["k_pb_off"], p["k_ps_on"], p["k_ps_off"]),
                                 (.005, .001, .01, .001))
                self.assertEqual((p["production_rate"], p["response_initial"],
                                  p["daughter_expression"], p["response_midpoint"]),
                                 (.005, 0, 0, 15))
                self.assertEqual((p["ps_crowding_threshold"], p["population_cap"]), (10, 80))

    def test_first_pair_is_an_exact_subset_of_standard_full_requests(self):
        pair = make_plan(json.loads((ROOT/"examples/bb170_standard_pair.json").read_text()))
        full = make_plan(json.loads((ROOT/"examples/bb170_standard_full.json").read_text()))
        self.assertEqual(len(pair), 2)
        self.assertEqual({j["seed"] for j in pair}, {2026})
        self.assertEqual({j["parameters"]["polymers"] for j in pair}, {0, 1000})
        by_id = {j["id"]: j["parameters"] for j in full}
        for job in pair:
            self.assertEqual(job["parameters"], by_id[job["id"]])

    def test_faster_growth_cli_records_boundary_and_labels_the_condition(self):
        # Very short, signal-free CLI fixture checks presentation and provenance.
        # It is not an assignment-scale BB170 simulation.
        with tempfile.TemporaryDirectory() as temporary:
            with contextlib.redirect_stdout(io.StringIO()):
                code = main(["--mode", "full", "--grid-size", "8", "--bacteria", "1",
                             "--strain", "BB170", "--duration", "0.01",
                             "--doubling-time", "10000", "--division-width", "20000",
                             "--production-rate", "0", "--d-s", "0", "--d-b", "0",
                             "--spatial-frames", "2", "--output", temporary])
            self.assertEqual(code, 0)
            run = next(Path(temporary).iterdir())
            metadata = json.loads((run/"run.json").read_text())
            self.assertEqual(metadata["population_model"]["division_support"], [0, 20000])
            self.assertIn("POP-07", metadata["interpretations"])
            drawing = (run/"population.svg").read_text()
            self.assertIn("BB170 faster growth: T=10000 a.u.", drawing)
            self.assertNotIn("Division timing differs from source", drawing)

    def test_faster_growth_ensemble_label_and_singleton_uncertainty(self):
        for seeds in [[21], [21, 22]]:
            with self.subTest(seeds=seeds), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                spec = small_spec()
                spec["seeds"] = seeds
                spec["base_parameters"].update(doubling_time=10000, division_width=20000)
                jobs = make_plan(spec)
                software = software_metadata()
                for job in jobs:
                    path = fixture_run(root, job, software)
                    job.update(status="complete", run_directory=str(path.relative_to(root)))
                save_ensemble(root, {"jobs": jobs, "software": software})
                drawing = (root/"ensemble.svg").read_text()
                self.assertIn("BB170 faster growth: T=10000 a.u.", drawing)
                self.assertNotIn("Division timing differs from source", drawing)
                if len(seeds) == 1:
                    self.assertIn("means only (n=1)", drawing)


if __name__ == "__main__":
    unittest.main()

