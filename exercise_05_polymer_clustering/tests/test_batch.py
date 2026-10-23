# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 20 October 2026

"""Check replicate identities, scientific aggregation, failures, and process isolation."""

import contextlib
import copy
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from polymer_clustering.batch import (
    available_cpus, controller_lock, execute_jobs, inspect_run, main, make_plan,
    matching_batch_runs, plan_hash, reuse_runs,
)
from polymer_clustering.batch_output import save_ensemble, statistics, write_csv
from polymer_clustering.output import software_metadata, write_json


def small_spec():
    """Accelerated BB170 is a software diagnostic, not an assignment result."""
    return {"base_parameters": {
        "grid_size": 8, "bacteria": 1, "strain": "BB170", "duration": 2,
        "sample_interval": 1, "spatial_frames": 2, "max_events": 10000,
        "doubling_time": 1, "division_width": 0, "population_cap": 4,
        "production_rate": 10, "k_pb_on": .005, "k_pb_off": .001,
        "k_ps_on": .01, "k_ps_off": .001, "font_family": "DejaVu Sans",
        "light_map_max": 5000,
    }, "conditions": [{"name": "control", "overrides": {"polymers": 0}},
                       {"name": "p3", "overrides": {"polymers": 4}}], "seeds": [21, 22]}


def fixture_run(root, job, software):
    """Synthetic, balanced observations make the aggregation arithmetic independent."""
    path = root / job["id"] / "fixture"
    path.mkdir(parents=True)
    p = job["parameters"]
    second_seed = job["seed"] == 22
    cells = 4 if second_seed else 2
    mean = 30 if second_seed else 10
    if job["condition"] == "p3":
        mean -= 8 if second_seed else 2
    response, binding, population = [], [], []
    for time in (0., 1., 2.):
        n = p["bacteria"] if time == 0 else cells
        s = 0 if time == 0 else n * 10
        light = 0 if time == 0 else mean
        produced = 0 if time == 0 else 100
        response.append(dict(time=time, bacteria=n, sensed_total=s, sensed_mean=s/n,
                             expression_total=s, expression_mean=s/n,
                             luminescence_total=n*light, luminescence_mean=light))
        binding.append(dict(time=time, free_signals=produced-s, bacteria_bound_signals=s,
                            polymer_bound_signals=0, total_signals=produced,
                            free_polymers=p["polymers"], attached_polymers=0,
                            pb_bonds=0, bridges=0, largest_cluster=1))
        population.append(dict(time=time, bacteria=n, divisions=n-p["bacteria"],
                               blocked_attempts=0, produced_signals=produced,
                               expected_total_signals=produced, actual_total_signals=produced))
    for name, rows in (("response", response), ("binding", binding), ("population", population)):
        write_csv(path / f"{name}_trajectory.csv", rows)
    write_json(path / "run.json", dict(status="complete", parameters=p, software=software,
               binding=dict(stop_time=2, events=100), response=response[-1],
               timing=dict(total_wall_seconds=10, cpu_seconds=9)))
    return path


class BatchChecks(unittest.TestCase):
    def test_unique_stable_seeds_and_input_validation(self):
        spec = small_spec()
        jobs = make_plan(spec)
        self.assertEqual([j["seed"] for j in jobs], [21, 22, 21, 22])
        self.assertEqual(len({j["id"] for j in jobs}), 4)
        changed = copy.deepcopy(jobs)
        changed[0]["status"] = "running"
        self.assertEqual(plan_hash(jobs), plan_hash(changed))
        changed[0]["parameters"]["seed"] = 23
        self.assertNotEqual(plan_hash(jobs), plan_hash(changed))
        for alter in (
            lambda s: s.update(seeds=[21, 21]),
            lambda s: s.update(seeds=[True]),
            lambda s: s["conditions"][1].update(name="../bad"),
            lambda s: s["conditions"][1].update(name="control"),
            lambda s: s["base_parameters"].update(k_ps_on=None),
            lambda s: s["base_parameters"].update(unknown=True),
            lambda s: s["conditions"][1]["overrides"].update(duration=3),
            lambda s: s["base_parameters"].update(mode="response"),
        ):
            invalid = copy.deepcopy(spec)
            alter(invalid)
            with self.assertRaises((ValueError, TypeError)):
                make_plan(invalid)

    def test_dry_run_and_invalid_workers_create_no_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = root / "batch.json"
            write_json(config, small_spec())
            output = root / "output"
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--config", str(config), "--output", str(output),
                                       "--dry-run"]), 0)
            self.assertFalse(output.exists())
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(["--config", str(config), "--output", str(output), "--workers", "0"])
            self.assertFalse(output.exists())

    def test_scheduler_cpu_limit(self):
        with patch.dict("os.environ", {"SLURM_CPUS_PER_TASK": "1"}):
            self.assertEqual(available_cpus(), 1)

    def test_only_one_controller_can_hold_a_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with controller_lock(root):
                with self.assertRaises(ValueError), controller_lock(root):
                    pass
            with controller_lock(root):
                pass

    def test_sample_sd_sem_and_singleton(self):
        result = statistics(np.array([1., 3.]))
        self.assertEqual(result["mean"], 2)
        self.assertAlmostEqual(result["sd"], np.sqrt(2))
        self.assertEqual(result["sem"], 1)
        self.assertIsNone(statistics(np.array([4.]))["sd"])
        self.assertIsNone(statistics(np.array([4.]))["sem"])

    def test_run_provenance_conservation_and_reuse(self):
        jobs, software = make_plan(small_spec()), software_metadata()
        with tempfile.TemporaryDirectory() as tmp:
            path = fixture_run(Path(tmp), jobs[0], software)
            meta, rows = inspect_run(path, jobs[0]["parameters"], software)
            self.assertEqual(len(rows), 3)
            # A previous release can be reused when its numerical sources agree.
            meta["software"]["program_version"] = "0.5.1"
            write_json(path / "run.json", meta)
            self.assertEqual(reuse_runs([path], jobs, software)[0][1]["id"], jobs[0]["id"])
            with self.assertRaises(ValueError):
                reuse_runs([path, path], jobs, software)
            meta["software"]["source_sha256"]["movement.py"] = "changed"
            write_json(path / "run.json", meta)
            with self.assertRaises(ValueError):
                inspect_run(path, jobs[0]["parameters"], software_metadata())
            meta["software"] = software_metadata()
            meta["status"] = "incomplete"
            write_json(path / "run.json", meta)
            with self.assertRaises(ValueError):
                inspect_run(path, jobs[0]["parameters"], software_metadata())
            meta["status"] = "complete"
            write_json(path / "run.json", meta)
            with (path / "binding_trajectory.csv").open() as f:
                binding = list(csv.DictReader(f))
            binding[-1]["total_signals"] = 99
            write_csv(path / "binding_trajectory.csv", binding)
            with self.assertRaises(ValueError):
                inspect_run(path, jobs[0]["parameters"], software_metadata())

    def test_equal_run_weights_paired_sd_and_no_partial_averages(self):
        jobs, software = make_plan(small_spec()), software_metadata()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for job in jobs:
                path = fixture_run(root, job, software)
                job.update(status="complete", run_directory=str(path.relative_to(root)))
            manifest = dict(jobs=jobs, software=software)
            save_ensemble(root, manifest)
            with (root / "ensemble_trajectory.csv").open() as f:
                rows = list(csv.DictReader(f))
            row = next(r for r in rows if r["condition"] == "control"
                       and r["time"] == "2.0" and r["metric"] == "luminescence_mean")
            # Mean(10,30)=20, not the cell-weighted (2*10+4*30)/6.
            self.assertEqual(float(row["mean"]), 20)
            self.assertAlmostEqual(float(row["sd"]), np.sqrt(200))
            with (root / "paired_differences.csv").open() as f:
                paired = list(csv.DictReader(f))
            row = next(r for r in paired if r["time"] == "2.0"
                       and r["metric"] == "luminescence_mean")
            self.assertEqual(float(row["mean"]), -5)
            self.assertAlmostEqual(float(row["sd"]), np.sqrt(18))
            jobs[0]["status"] = "incomplete"
            with self.assertRaises(ValueError):
                save_ensemble(root, manifest)

    def test_interrupt_terminates_children_and_preserves_attempt(self):
        jobs, software = make_plan(small_spec())[:1], software_metadata()
        class Child:
            returncode = None
            def poll(self): return None
            def terminate(self): self.returncode = -15
            def wait(self, timeout=None): return self.returncode
        child = Child()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = dict(jobs=jobs, software=software)
            with patch("polymer_clustering.batch.subprocess.Popen", return_value=child), \
                 patch("polymer_clustering.batch.time.sleep", side_effect=KeyboardInterrupt), \
                 contextlib.redirect_stdout(io.StringIO()), self.assertRaises(KeyboardInterrupt):
                execute_jobs(root, manifest, 1)
            self.assertEqual(child.returncode, -15)
            self.assertEqual(jobs[0]["status"], "interrupted")
            self.assertEqual(jobs[0]["attempts"][0]["exit_code"], -15)
            self.assertTrue((root / "batch.json").exists())

    def test_incomplete_batch_is_not_averaged_and_resume_retains_attempt(self):
        spec = small_spec()
        spec["base_parameters"]["max_events"] = 1
        spec["conditions"] = spec["conditions"][:1]
        spec["seeds"] = [21]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = root / "request.json"
            write_json(config, spec)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--config", str(config), "--output", str(root/"out")]), 3)
            directory = next((root/"out").iterdir())
            self.assertFalse((directory / "ensemble_trajectory.csv").exists())
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--resume", str(directory)]), 3)
            manifest = json.loads((directory/"batch.json").read_text())
            self.assertEqual(manifest["status"], "incomplete")
            self.assertEqual(len(manifest["jobs"][0]["attempts"]), 2)
            self.assertEqual(manifest["jobs"][0]["status"], "incomplete")
            self.assertTrue((directory/manifest["jobs"][0]["attempts"][0]["directory"]).exists())
            self.assertFalse((directory / "ensemble.png").exists())

    def test_imported_jobs_are_copied_and_skip_execution(self):
        spec = small_spec()
        spec["seeds"] = [21]
        spec["conditions"] = spec["conditions"][:1]
        jobs, software = make_plan(spec), software_metadata()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = fixture_run(root/"original", jobs[0], software)
            config = root / "request.json"
            # A changed completed-run event budget is harmless and stays recorded.
            spec["base_parameters"]["max_events"] = 20000
            write_json(config, spec)
            with patch("polymer_clustering.batch.subprocess.Popen") as launch, \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--config", str(config), "--reuse-run", str(source),
                                       "--output", str(root/"out")]), 0)
                launch.assert_not_called()
            directory = next((root/"out").iterdir())
            manifest = json.loads((directory/"batch.json").read_text())
            job = manifest["jobs"][0]
            self.assertEqual(job["reused_from"], str(source))
            self.assertEqual(job["attempts"], [])
            self.assertTrue((source/"run.json").exists())
            self.assertEqual((directory/job["run_directory"]/"run.json").read_bytes(),
                             (source/"run.json").read_bytes())

    def test_reuse_batch_matches_settings_preserves_originals_and_detects_duplicates(self):
        spec = small_spec()
        jobs, software = make_plan(spec), software_metadata()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = root/"original_batch"
            for job in jobs:
                run = fixture_run(original, job, software)
                job.update(status="complete", run_directory=str(run.relative_to(original)))
            old = dict(schema_version=1, status="complete", jobs=jobs,
                       plan_sha256=plan_hash(jobs), software=software)
            write_json(original/"batch.json", old)
            # A new label does not alter a numerical request. An extra condition
            # remains pending, rather than being accidentally matched by seed.
            expanded = copy.deepcopy(spec)
            expanded["conditions"][1]["name"] = "renamed_p3"
            expanded["conditions"].append({"name": "pb_only", "overrides": {
                "polymers": 4, "k_ps_on": 0,
            }})
            target = make_plan(expanded)
            selected = matching_batch_runs([original], target)
            imports = reuse_runs(selected, target, software)
            self.assertEqual(len(imports), 4)
            self.assertEqual({j["condition"] for _, j in imports}, {"control", "renamed_p3"})
            with self.assertRaises(ValueError):
                reuse_runs(selected+selected, target, software)
            # The CLI also supports importing an entire compatible complete batch.
            config = root/"request.json"
            write_json(config, spec)
            with patch("polymer_clustering.batch.subprocess.Popen") as launch, \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--config", str(config), "--reuse-batch", str(original),
                                       "--output", str(root/"new")]), 0)
                launch.assert_not_called()
            new = next((root/"new").iterdir())
            saved = json.loads((new/"batch.json").read_text())
            self.assertEqual(saved["reuse_source_batches"], [str(original)])
            self.assertTrue(all(j["status"] == "complete" and not j["attempts"]
                                for j in saved["jobs"]))
            self.assertEqual(json.loads((original/"batch.json").read_text()), old)

    def test_reuse_batch_rejects_incomplete_mismatched_and_escaping_records(self):
        spec = small_spec()
        jobs, software = make_plan(spec), software_metadata()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for job in jobs:
                run = fixture_run(root, job, software)
                job.update(status="complete", run_directory=str(run.relative_to(root)))
            previous = dict(schema_version=1, status="incomplete", jobs=jobs,
                            plan_sha256=plan_hash(jobs))
            write_json(root/"batch.json", previous)
            with self.assertRaises(ValueError):
                matching_batch_runs([root], make_plan(spec))
            previous["status"] = "complete"
            write_json(root/"batch.json", previous)
            different = copy.deepcopy(spec)
            different["base_parameters"]["duration"] = 3
            with self.assertRaises(ValueError):
                matching_batch_runs([root], make_plan(different))
            previous["jobs"][0]["run_directory"] = "../outside"
            write_json(root/"batch.json", previous)
            with self.assertRaises(ValueError):
                matching_batch_runs([root], make_plan(spec))

    def test_pairwise_contrasts_follow_seed_pairing_and_direction(self):
        spec = small_spec()
        spec["conditions"].append({"name": "pb_only", "overrides": {
            "polymers": 4, "k_ps_on": 0,
        }})
        jobs, software = make_plan(spec), software_metadata()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for job in jobs:
                run = fixture_run(root, job, software)
                job.update(status="complete", run_directory=str(run.relative_to(root)))
            save_ensemble(root, dict(jobs=jobs, software=software))
            with (root/"pairwise_differences.csv").open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 3*3*11)
            row = next(r for r in rows if r["condition"] == "pb_only"
                       and r["reference"] == "p3" and r["time"] == "2.0"
                       and r["metric"] == "luminescence_mean")
            # Synthetic p3 is lower by 2 and 8, so PB-only minus P3 is [2,8].
            self.assertEqual(float(row["mean"]), 5)
            self.assertAlmostEqual(float(row["sd"]), np.sqrt(18))

    def test_real_serial_parallel_replay_and_resume(self):
        # Actual separate interpreters check that scheduling never selects seeds.
        package_root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = root / "request.json"
            write_json(config, small_spec())
            batches = []
            for workers in (1, 2) if available_cpus() >= 2 else (1,):
                output = root / str(workers)
                process = subprocess.run(
                    [sys.executable, "-m", "polymer_clustering.batch", "--config", str(config),
                     "--output", str(output), "--workers", str(workers)], cwd=package_root,
                    capture_output=True, text=True, timeout=180,
                )
                self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
                batches.append(next(output.iterdir()))
            first = json.loads((batches[0] / "batch.json").read_text())
            if len(batches) == 2:
                second = json.loads((batches[1] / "batch.json").read_text())
                for a, b in zip(first["jobs"], second["jobs"]):
                    pa, pb = batches[0] / a["run_directory"], batches[1] / b["run_directory"]
                    for name in pa.glob("*.csv"):
                        self.assertEqual(name.read_bytes(), (pb / name.name).read_bytes(), name.name)
                    for name in pa.glob("*.npz"):
                        with np.load(name) as aa, np.load(pb/name.name) as bb:
                            self.assertEqual(aa.files, bb.files)
                            for key in aa.files:
                                np.testing.assert_array_equal(aa[key], bb[key], err_msg=name.name)
                    ma, mb = [json.loads((p / "run.json").read_text()) for p in (pa, pb)]
                    for key in ("response", "population", "binding", "movement_event_counts"):
                        self.assertEqual(ma[key], mb[key], key)
                for name in ("ensemble_trajectory.csv", "paired_differences.csv"):
                    self.assertEqual((batches[0] / name).read_bytes(), (batches[1] / name).read_bytes())
            # Complete jobs are inspected and retained without launching children.
            with patch("polymer_clustering.batch.subprocess.Popen") as launch, \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["--resume", str(batches[0])]), 0)
                launch.assert_not_called()
            final = json.loads((batches[0] / "batch.json").read_text())
            self.assertEqual([len(j["attempts"]) for j in final["jobs"]], [1]*4)


if __name__ == "__main__":
    unittest.main()
