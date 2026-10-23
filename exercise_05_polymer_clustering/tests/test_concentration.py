# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 21 October 2026

"""Keep the focused concentration pilot matched to the completed affinity plan."""

import json
from pathlib import Path
import unittest

from polymer_clustering.batch import comparable_parameters, make_plan


ROOT = Path(__file__).resolve().parents[1]


class ConcentrationChecks(unittest.TestCase):
    def test_eight_jobs_reuse_four_original_requests_without_presentation_drift(self):
        pilot = make_plan(json.loads((ROOT/"examples/mm32_concentration_pilot.json").read_text()))
        affinity = make_plan(json.loads((ROOT/"examples/mm32_affinity_full.json").read_text()))
        self.assertEqual(len(pilot), 8)
        self.assertEqual({j["seed"] for j in pilot}, {2026, 2027})
        reused = [job for job in pilot if any(
            comparable_parameters(job["parameters"]) == comparable_parameters(old["parameters"])
            for old in affinity)]
        self.assertEqual({j["id"] for j in reused}, {
            f"{name}/seed_{seed}" for name in ("no_polymer", "p3_1000") for seed in (2026, 2027)})
        self.assertEqual(len(pilot) - len(reused), 4)

    def test_polymer_count_is_the_only_changed_input_within_each_seed(self):
        jobs = make_plan(json.loads((ROOT/"examples/mm32_concentration_pilot.json").read_text()))
        for seed in (2026, 2027):
            parameters = [j["parameters"] for j in jobs if j["seed"] == seed]
            self.assertEqual([p["polymers"] for p in parameters], [0, 1000, 3000, 7000])
            for setting in parameters:
                self.assertEqual({k: v for k, v in setting.items() if k != "polymers"},
                                 {k: v for k, v in parameters[0].items() if k != "polymers"})
                self.assertEqual(setting["duration"], 30000)
                self.assertEqual(setting["ps_crowding_threshold"], 10)
                self.assertEqual(setting["spatial_frames"], 61)
                self.assertEqual(setting["response_initial"], 0)
                self.assertEqual(setting["font_family"], "Arial")


if __name__ == "__main__":
    unittest.main()
