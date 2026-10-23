# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 21 October 2026

"""Check spatial arithmetic against known states and preserve completed input runs."""

import contextlib
import copy
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

import numpy as np

from polymer_clustering.batch import make_plan, plan_hash
from polymer_clustering.batch_output import write_csv
from polymer_clustering.binding import BindingConfig, BindingState, polymer_reactions, sample_counts
from polymer_clustering.model import Bacterium
from polymer_clustering.output import software_metadata, write_json
from polymer_clustering.spatial import (
    aggregate, analyse, file_hash, load_state, main, measure_state, select_jobs,
)


def known_state():
    """Ten polymers crowd one cell site; one other polymer remains outside cells.

    Three signals are polymer-bound: two at the crowded site and one outside.
    Thus 48 empty PS slots are blocked and four remain accessible. Two free
    signals sit at the crowded site; the other 9995 are outside the cell.
    """
    occupied = np.full((4, 4), -1, dtype=np.int32)
    occupied[:2, :2] = 0
    free = np.zeros((4, 4), dtype=np.int64)
    free[0, 0], free[3, 3] = 2, 9995
    state = BindingState(
        occupancy=occupied, bacteria=(Bacterium(0, 0, 0),), free_signals=free,
        bacteria_signals=np.zeros((4, 4), dtype=np.int64),
        bacteria_polymers=np.zeros((4, 4), dtype=np.int64),
        polymer_sites=np.array([0] * 10 + [15], dtype=np.int64),
        polymer_signals=np.array([2] + [0] * 9 + [1], dtype=np.int64),
        polymer_contacts=np.full((11, 2), -1, dtype=np.int64),
    )
    state.validate(10000, 11)
    return state


def save_arrays(path, state):
    """Use the same numeric-only archive schema as the real binding output."""
    arrays = {k: v for k, v in vars(state).items() if k != "bacteria"}
    arrays["bacteria_positions"] = np.array([(b.identifier, b.x, b.y) for b in state.bacteria])
    np.savez_compressed(path, **arrays)


def fixture_batch(root):
    """Create balanced synthetic histories to exercise the read-only CLI boundary.

    These observations are arithmetic fixtures, not scientifically simulated
    trajectories. Actual saved scientific runs are checked separately.
    """
    jobs = make_plan({"base_parameters": {
        "grid_size": 4, "bacteria": 1, "polymers": 11, "duration": 1,
        "sample_interval": 1, "k_pb_on": .005, "k_pb_off": .001,
        "k_ps_on": .01, "k_ps_off": .001, "font_family": "DejaVu Sans",
    }, "conditions": [{"name": "without_pb", "overrides": {"k_pb_on": 0}},
                       {"name": "with_pb", "overrides": {}}], "seeds": [4, 5]})
    software = software_metadata()
    for job in jobs:
        directory = root / job["id"] / "saved_run"
        directory.mkdir(parents=True)
        state = known_state()
        if job["condition"] == "without_pb":
            state.polymer_sites = np.arange(11, dtype=np.int64)
        save_arrays(directory / "binding_state.npz", state)
        binding = [sample_counts(state, t) for t in (0., 1.)]
        light = 10 * job["seed"] + (2 if job["condition"] == "with_pb" else 0)
        response = [dict(time=t, bacteria=1, sensed_total=0, sensed_mean=0,
                         expression_total=0, expression_mean=0,
                         luminescence_total=light, luminescence_mean=light) for t in (0., 1.)]
        population = [dict(time=t, bacteria=1, divisions=0, blocked_attempts=0,
                           produced_signals=0, expected_total_signals=10000,
                           actual_total_signals=10000) for t in (0., 1.)]
        for name, rows in (("binding", binding), ("response", response), ("population", population)):
            write_csv(directory / f"{name}_trajectory.csv", rows)
        old = copy.deepcopy(software)
        old["program_version"] = "0.5.1"
        write_json(directory / "run.json", dict(status="complete", parameters=job["parameters"],
                   software=old, binding=dict(stop_time=1), response=response[-1]))
        job.update(status="complete", run_directory=str(directory.relative_to(root)))
    manifest = dict(schema_version=1, status="complete", jobs=jobs,
                    plan_sha256=plan_hash(jobs), software=software)
    write_json(root / "batch.json", manifest)
    return manifest


class SpatialChecks(unittest.TestCase):
    def test_known_counts_capacity_and_instantaneous_rates(self):
        state = known_state()
        values, maps = measure_state(state, {"k_ps_on": .01, "k_ps_off": .001})
        self.assertEqual(values["crowded_sites"], 1)
        self.assertEqual(values["polymers_at_crowded_sites"], 10)
        self.assertEqual(values["polymers_on_bacterial_sites"], 10)
        self.assertEqual(values["empty_ps_slots"], 52)
        self.assertEqual(values["blocked_empty_ps_slots"], 48)
        self.assertAlmostEqual(values["ps_capacity_occupancy"], 3 / 55)
        self.assertAlmostEqual(values["ps_bind_propensity"], 399.8)
        self.assertAlmostEqual(values["ps_bind_propensity_without_cutoff"], 400.76)
        self.assertAlmostEqual(values["ps_bind_propensity_lost_to_cutoff"], .96)
        self.assertAlmostEqual(values["ps_unbind_propensity"], .003)
        self.assertEqual(int(maps["polymer_bound_signals"].sum()), 3)
        # Independent check: the actual reaction-channel construction agrees
        # with this endpoint diagnostic's predicted aggregate PS on rate.
        config = BindingConfig(k_pb_on=0, k_pb_off=.001, k_ps_on=.01, k_ps_off=.001)
        actual = sum(channel.propensity for p in range(11)
                     for channel in polymer_reactions(state, config, {}, p,
                         int(maps["polymer_density"].ravel()[state.polymer_sites[p]]))
                     if channel.kind == "ps_bind")
        self.assertAlmostEqual(actual, values["ps_bind_propensity"])

    def test_cutoff_boundary_and_previously_bound_signals(self):
        state = known_state()
        state.polymer_sites[1] = 1  # Nine at the original anchor no longer meet P >= 10.
        values, _ = measure_state(state, {"k_ps_on": .01, "k_ps_off": .001})
        self.assertEqual(values["crowded_sites"], 0)
        self.assertEqual(values["blocked_empty_ps_slots"], 0)
        self.assertEqual(values["ps_bind_propensity_lost_to_cutoff"], 0)
        self.assertAlmostEqual(values["ps_unbind_propensity"], .003)
        self.assertEqual(int(state.polymer_signals.sum()), 3)

    def test_zero_polymer_ratios_and_zero_affinity(self):
        state = known_state()
        state.polymer_sites = state.polymer_signals = np.empty(0, dtype=np.int64)
        state.polymer_contacts = np.empty((0, 2), dtype=np.int64)
        values, _ = measure_state(state, {})
        for key in ("ps_capacity_occupancy", "polymer_crowded_fraction",
                    "polymer_on_bacteria_fraction", "blocked_empty_ps_fraction"):
            self.assertIsNone(values[key])
        self.assertEqual(values["blocked_empty_ps_slots"], 0)
        values, _ = measure_state(known_state(), {"k_ps_on": 0, "k_ps_off": .001})
        self.assertEqual(values["ps_bind_propensity_without_cutoff"], 0)
        self.assertAlmostEqual(values["ps_unbind_propensity"], .003)

    def test_paired_statistics_use_seed_differences_and_missing_ratios(self):
        rows = [{"condition": c, "seed": seed, "run_directory": "unused", "time": 1,
                 "value": value, "ratio": None}
                for c, seed, value in (("a", 4, 10), ("a", 5, 30), ("b", 4, 8), ("b", 5, 22))]
        ensemble, paired = aggregate(rows)
        result = next(r for r in paired if r["metric"] == "value")
        self.assertEqual(result["mean"], -5)
        self.assertAlmostEqual(result["sd"], np.sqrt(18))
        self.assertEqual(next(r for r in ensemble if r["metric"] == "ratio")["n"], 0)
        self.assertIsNone(aggregate(rows[:1])[0][0]["sd"])

    def test_final_array_count_and_integer_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "state.npz"
            state = known_state()
            save_arrays(path, state)
            final = {**sample_counts(state, 1), "bacteria": 1}
            parameters = {"grid_size": 4, "polymers": 11}
            load_state(path, parameters, final)
            with self.assertRaises(ValueError):
                load_state(path, parameters, {**final, "polymer_bound_signals": 4})
            state.polymer_signals = state.polymer_signals.astype(float)
            save_arrays(path, state)
            with self.assertRaises(ValueError):
                load_state(path, parameters, final)

    def test_plan_completion_seed_balance_and_unknown_selections(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = fixture_batch(Path(tmp))
            self.assertEqual(len(select_jobs(manifest, None, None)), 4)
            self.assertEqual(len(select_jobs(manifest, ["with_pb"], [4])), 1)
            for names, seeds in ((["unknown"], None), (["with_pb"], [99]),
                                 (["with_pb", "with_pb"], None)):
                with self.assertRaises(ValueError):
                    select_jobs(manifest, names, seeds)
            bad = copy.deepcopy(manifest)
            bad["jobs"][0]["status"] = "incomplete"
            with self.assertRaises(ValueError):
                select_jobs(bad, None, None)
            bad = copy.deepcopy(manifest)
            bad["jobs"][0]["parameters"]["k_ps_on"] = .1
            with self.assertRaises(ValueError):
                select_jobs(bad, None, None)
            with self.assertRaises(ValueError):
                analyse(Path(tmp), None, None, 99)

    def test_outside_paths_rejected_before_writing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = fixture_batch(root)
            manifest["jobs"][0]["run_directory"] = "../outside"
            write_json(root / "batch.json", manifest)
            with self.assertRaises(ValueError):
                analyse(root, None, None, None)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--batch", str(root), "--output", str(root / "analysis")]), 1)
            self.assertFalse((root / "analysis").exists())

    def test_cli_preserves_inputs_and_exports_compact_bundle_without_simulating(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, output = Path(tmp) / "batch", Path(tmp) / "output"
            fixture_batch(root)
            before = {p: file_hash(p) for p in root.rglob("*") if p.is_file()}
            args = ["--batch", str(root), "--output", str(output), "--font-family", "DejaVu Sans"]
            with contextlib.redirect_stdout(io.StringIO()), \
                 patch("polymer_clustering.movement.simulate_movement", side_effect=AssertionError):
                self.assertEqual(main([*args, "--dry-run"]), 0)
                self.assertFalse(output.exists())
                self.assertEqual(main(args), 0)
            self.assertEqual(before, {p: file_hash(p) for p in root.rglob("*") if p.is_file()})
            directory = next(output.iterdir())
            report = json.loads((directory / "spatial_report.json").read_text())
            self.assertEqual(report["status"], "complete")
            self.assertEqual(report["map_seed"], 4)
            with (directory / "spatial_sites.csv").open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), 64)
            with ZipFile(directory / "spatial_summary.zip") as archive:
                self.assertEqual(len(archive.namelist()), 6)
                self.assertNotIn("spatial_sites.csv", archive.namelist())
                self.assertIsNone(archive.testzip())


if __name__ == "__main__":
    unittest.main()
