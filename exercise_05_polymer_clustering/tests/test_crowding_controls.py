# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 21 October 2026

"""Verify the PS-only intervention, old-default compatibility and child-process settings."""

import contextlib
import copy
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from polymer_clustering.batch import (
    LEGACY_BINDING_SHA256, PARAMETERISED_BINDING_SHA256, available_cpus,
    inspect_run, main as batch_main, make_plan, numerical_sources_match,
)
from polymer_clustering.binding import BindingConfig, polymer_reactions, reaction_channels
from polymer_clustering.cli import build_parser, main as run_main
from polymer_clustering.movement import MovementConfig, MovementEngine
from polymer_clustering.output import software_metadata, write_json
from polymer_clustering.spatial import measure_state
from test_spatial import fixture_batch, known_state


ROOT = Path(__file__).resolve().parents[1]


class CrowdingChecks(unittest.TestCase):
    def test_threshold_validation_and_explicit_off_flag(self):
        self.assertEqual(BindingConfig().ps_crowding_threshold, 10)
        self.assertIsNone(BindingConfig(ps_crowding_threshold=None).ps_crowding_threshold)
        for value in (0, -1, True, 10., "off"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                BindingConfig(ps_crowding_threshold=value)
        parser = build_parser()
        self.assertFalse(hasattr(parser.parse_args([]), "ps_crowding_threshold"))
        self.assertIsNone(parser.parse_args(["--ps-crowding-threshold", "off"]).ps_crowding_threshold)
        self.assertEqual(parser.parse_args(["--ps-crowding-threshold", "20"]).ps_crowding_threshold, 20)

    def test_cutoff_changes_ps_on_only_and_retains_release(self):
        state = known_state()
        common = dict(k_pb_on=.005, k_pb_off=.001, k_ps_on=.01, k_ps_off=.001)
        standard = BindingConfig(**common)
        relaxed = BindingConfig(**common, ps_crowding_threshold=20)
        disabled = BindingConfig(**common, ps_crowding_threshold=None)
        channels = [polymer_reactions(state, c, {}, 0, 10) for c in (standard, relaxed, disabled)]
        self.assertFalse(any(c.kind == "ps_bind" for c in channels[0]))
        self.assertEqual(next(c.propensity for c in channels[1] if c.kind == "ps_bind"), .06)
        self.assertEqual(channels[1], channels[2])
        self.assertEqual([c for c in channels[0] if c.kind != "ps_bind"],
                         [c for c in channels[2] if c.kind != "ps_bind"])
        self.assertEqual(next(c.propensity for c in channels[0] if c.kind == "ps_unbind"), .002)
        self.assertFalse(any(c.kind == "ps_bind" for c in polymer_reactions(state, relaxed, {}, 0, 20)))
        self.assertTrue(any(c.kind == "ps_bind" for c in polymer_reactions(state, disabled, {}, 0, 100)))

    def test_cached_and_reference_channels_agree_for_each_cutoff(self):
        for threshold in (10, 20, None):
            with self.subTest(threshold=threshold):
                state = known_state()
                config = BindingConfig(k_pb_on=.005, k_pb_off=.001, k_ps_on=.01,
                                       k_ps_off=.001, ps_crowding_threshold=threshold)
                engine = MovementEngine(state, config, MovementConfig())
                expected = reaction_channels(state, config, engine.neighbours)
                self.assertAlmostEqual(engine.bs_tree.total + engine.chem_tree.total,
                                       sum(c.propensity for c in expected))

    def test_spatial_metrics_distinguish_density_from_disabled_occlusion(self):
        values, _ = measure_state(known_state(), {"k_ps_on": .01, "k_ps_off": .001,
                                                   "ps_crowding_threshold": None})
        self.assertEqual(values["crowded_sites"], 0)
        self.assertEqual(values["blocked_empty_ps_slots"], 0)
        self.assertEqual(values["sites_at_source_density_threshold"], 1)
        self.assertAlmostEqual(values["polymer_dense_fraction_at_10"], 10 / 11)
        self.assertAlmostEqual(values["ps_bind_propensity"], 400.76)
        self.assertEqual(values["ps_bind_propensity_lost_to_cutoff"], 0)

    def test_pilot_is_a_matched_single_rule_intervention(self):
        jobs = make_plan(json.loads((ROOT / "examples/mm32_crowding_pilot.json").read_text()))
        self.assertEqual(len(jobs), 4)
        for seed in (2026, 2027):
            pair = [j["parameters"] for j in jobs if j["seed"] == seed]
            self.assertEqual({k for k in pair[0] if pair[0][k] != pair[1][k]},
                             {"ps_crowding_threshold"})
            self.assertEqual(pair[0]["ps_crowding_threshold"], 10)
            self.assertIsNone(pair[1]["ps_crowding_threshold"])
            self.assertEqual(pair[0]["duration"], 10000)
            self.assertEqual(pair[0]["k_pb_on"], .005)
            self.assertEqual(pair[0]["k_ps_on"], .01)
            self.assertTrue(pair[0]["growth"])

    def test_legacy_bridge_is_limited_to_exact_hashes_and_threshold_ten(self):
        current = software_metadata()
        self.assertEqual(current["source_sha256"]["binding.py"], PARAMETERISED_BINDING_SHA256)
        old = copy.deepcopy(current)
        old["source_sha256"]["binding.py"] = LEGACY_BINDING_SHA256
        self.assertTrue(numerical_sources_match(old, current, {}))
        for threshold in (None, 20):
            self.assertFalse(numerical_sources_match(old, current, {"ps_crowding_threshold": threshold}))
        edited = copy.deepcopy(current)
        edited["source_sha256"]["binding.py"] = "unreviewed edit"
        self.assertFalse(numerical_sources_match(old, edited, {}))
        changed = copy.deepcopy(old)
        changed["source_sha256"]["movement.py"] = "different transport"
        self.assertFalse(numerical_sources_match(changed, current, {}))
        changed = copy.deepcopy(old)
        changed["numpy"] = "different NumPy"
        self.assertFalse(numerical_sources_match(changed, current, {}))

    def test_old_missing_setting_is_default_ten_not_disabled(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = fixture_batch(root)
            job = manifest["jobs"][0]
            path = root / job["run_directory"]
            meta = json.loads((path / "run.json").read_text())
            meta["parameters"].pop("ps_crowding_threshold")
            meta["software"]["source_sha256"]["binding.py"] = LEGACY_BINDING_SHA256
            write_json(path / "run.json", meta)
            inspect_run(path, job["parameters"], software_metadata())
            with self.assertRaises(ValueError):
                inspect_run(path, {**job["parameters"], "ps_crowding_threshold": None}, software_metadata())

    def test_cli_explicit_off_overrides_json_and_records_intervention(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = root / "config.json"
            write_json(config, {"mode": "binding", "grid_size": 4, "bacteria": 1,
                               "polymers": 11, "seed": 17, "rng_scheme": "independent_v1",
                               "duration": .1, "sample_interval": .1, "k_pb_on": .005,
                               "k_pb_off": .001, "k_ps_on": .01, "k_ps_off": .001,
                               "ps_crowding_threshold": 20, "font_family": "DejaVu Sans"})
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(run_main(["--config", str(config), "--ps-crowding-threshold", "off",
                                           "--output", str(root / "out")]), 0)
            meta = json.loads(next((root / "out").glob("*/run.json")).read_text())
            self.assertIsNone(meta["parameters"]["ps_crowding_threshold"])
            self.assertFalse(meta["ps_crowding_rule"]["matches_source_default"])
            self.assertIn("CROWD-01", meta["interpretations"])

    def test_null_setting_survives_actual_batch_child_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            spec = json.loads((ROOT / "examples/mm32_crowding_pilot.json").read_text())
            spec["base_parameters"].update(grid_size=6, bacteria=2, polymers=360, duration=2,
                                           sample_interval=1, spatial_frames=2,
                                           font_family="DejaVu Sans")
            spec["seeds"] = [17]
            write_json(root / "config.json", spec)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(batch_main(["--config", str(root / "config.json"), "--workers",
                                            str(min(2, available_cpus())), "--output", str(root / "out")]), 0)
            batch = next((root / "out").iterdir())
            manifest = json.loads((batch / "batch.json").read_text())
            for job in manifest["jobs"]:
                meta = json.loads((batch / job["run_directory"] / "run.json").read_text())
                self.assertEqual(meta["parameters"]["ps_crowding_threshold"],
                                 job["parameters"]["ps_crowding_threshold"])
                self.assertEqual(meta["status"], "complete")
            before = [j["attempts"] for j in manifest["jobs"]]
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(batch_main(["--resume", str(batch)]), 0)
            after = json.loads((batch / "batch.json").read_text())
            self.assertEqual(before, [j["attempts"] for j in after["jobs"]])


if __name__ == "__main__":
    unittest.main()
