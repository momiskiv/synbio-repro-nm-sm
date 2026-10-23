# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 20 October 2026

"""Check source affinity ratios, matched experimental settings and disabled reactions."""

import json
from pathlib import Path
import unittest

from polymer_clustering.batch import make_plan
from polymer_clustering.binding import BindingConfig
from polymer_clustering.model import InitialConfig, initialise
from polymer_clustering.movement import MovementConfig, simulate_movement
from polymer_clustering.population import PopulationConfig
from polymer_clustering.response import ResponseConfig


ROOT = Path(__file__).resolve().parents[1]


class AffinityChecks(unittest.TestCase):
    def test_source_ratios_and_explicit_zero_forward_rates(self):
        spec = json.loads((ROOT/"examples/mm32_affinity_full.json").read_text())
        expected = {"no_polymer": (5, 10), "p1_like_1000": (5, 0),
                    "p2_like_1000": (0, .1), "p3_no_pb_1000": (0, 10),
                    "p3_1000": (5, 10)}
        self.assertEqual(spec["seeds"], list(range(2026, 2036)))
        jobs = make_plan(spec)
        self.assertEqual(len(jobs), 50)
        base = jobs[0]["parameters"]
        for job in jobs:
            p = job["parameters"]
            actual = (p["k_pb_on"]/p["k_pb_off"], p["k_ps_on"]/p["k_ps_off"])
            for value, target in zip(actual, expected[job["condition"]]):
                self.assertAlmostEqual(value, target)
            self.assertEqual(p["k_pb_off"], .001)
            self.assertEqual(p["k_ps_off"], .001)
            self.assertEqual(p["polymers"], 0 if job["condition"] == "no_polymer" else 1000)
            # Affinity controls must not accidentally change transport, growth,
            # receptor response, sampling, event budget, or initial population.
            for key in p.keys() - {"seed", "polymers", "k_pb_on", "k_ps_on"}:
                self.assertEqual(p[key], base[key], key)

    def test_original_comparison_remains_exactly_reusable(self):
        original = make_plan(json.loads((ROOT/"examples/mm32_replicates.json").read_text()))
        extended = make_plan(json.loads((ROOT/"examples/mm32_affinity_full.json").read_text()))
        lookup = {j["id"]: j for j in extended}
        for job in original:
            self.assertEqual(job["parameters"], lookup[job["id"]]["parameters"])
        pilot = make_plan(json.loads((ROOT/"examples/mm32_affinity_pilot.json").read_text()))
        self.assertEqual(len(pilot), 10)
        self.assertEqual({j["seed"] for j in pilot}, {2026, 2027})
        for job in pilot:
            full = lookup[job["id"]]["parameters"]
            for key in full.keys() - {"duration", "spatial_frames", "max_events"}:
                self.assertEqual(job["parameters"][key], full[key], key)
            self.assertEqual(job["parameters"]["duration"], 1000)

    def test_disabled_binding_classes_remain_absent_through_actual_dynamics(self):
        # A deliberately smaller grid makes these structural checks inexpensive.
        # These runs verify disabled reactions, not the long MM32 response curve.
        for name, pb_on, ps_on in (("p1", .005, 0), ("p2", 0, .0001),
                                   ("p3_without_pb", 0, .01)):
            with self.subTest(condition=name):
                state = initialise(InitialConfig(grid_size=8, bacteria=2, polymers=50, seed=17))
                result = simulate_movement(
                    state, BindingConfig(duration=10, sample_interval=5, max_events=100000,
                        k_pb_on=pb_on, k_pb_off=.001, k_ps_on=ps_on, k_ps_off=.001),
                    MovementConfig(spatial_frames=2), seed=17,
                    response=ResponseConfig(), strain="MM32", population=PopulationConfig(),
                )
                self.assertEqual(result.status, "complete")
                result.state.validate(10000, 50)
                if ps_on == 0:
                    self.assertTrue(all(r["polymer_bound_signals"] == 0 for r in result.trajectory))
                    self.assertEqual(result.event_counts.get("ps_bind", 0), 0)
                    self.assertGreater(result.event_counts.get("pb_bind", 0), 0)
                if pb_on == 0:
                    self.assertTrue(all(r["pb_bonds"] == r["bridges"] == 0
                                        and r["largest_cluster"] == 1 for r in result.trajectory))
                    self.assertEqual(result.event_counts.get("pb_bind", 0), 0)
                    self.assertGreater(result.event_counts.get("ps_bind", 0), 0)


if __name__ == "__main__":
    unittest.main()
