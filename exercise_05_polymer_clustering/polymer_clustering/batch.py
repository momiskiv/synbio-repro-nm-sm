# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Run independent, explicitly seeded simulations in isolated Python processes.

The controller only schedules whole runs. It never partitions an SSA trajectory
or draws random numbers. Completed jobs survive interruption; unfinished jobs
are replayed from their original seed when the batch is resumed.
"""

import argparse
from contextlib import contextmanager
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from .binding import BindingConfig, CROWDING_THRESHOLD
from .model import InitialConfig
from .movement import MovementConfig
from .output import software_metadata, write_json
from .plotting import PlotConfig, resolve_font
from .population import PopulationConfig
from .progress import BatchProgress, add_progress_argument
from .response import ResponseConfig


# These modules define the numerical reconstruction. Reuse requires matching
# NumPy and all six fingerprints, apart from the exact reviewed bridges below.
NUMERICAL_MODULES = ("model.py", "binding.py", "movement.py", "rates.py",
                     "response.py", "population.py")
# Reviewed progress and optimisation changes preserve numerical trajectories.
# Accept only these exact predecessor/target bytes, never release-name guesses.
# Legacy binding remains restricted to cutoff ten; all other modules and NumPy
# must match, and completed-run balances/grids are still independently checked.
LEGACY_BINDING_SHA256 = "22e480f462e9dab69033e1f465125ece08320952dfbddd5509cf602b0da91d94"
PRE_PROGRESS_BINDING_SHA256 = "17f7a817f07f8d395887762f6ed93e3fc735e86d571924dbe1f9fef0adbce807"
PARAMETERISED_BINDING_SHA256 = "1c5f02f0af61df5f05c10314fdcea5c7f2b767f7d2fca10b45c63040dcc7443f"
PRE_PROGRESS_MOVEMENT_SHA256 = "2de7680f1aeec7ba2da7a7bd706bd2f82d6f81efc21352362ce11b547f30fe66"
PROGRESS_MOVEMENT_SHA256 = "c563f5d3ad073d7c9791b2092f1ae9874f0f5072707c701685f27da5970e467d"
OPTIMISED_MOVEMENT_SHA256 = "7c49b22443bb309b1efbd8fdd0053b2e7f2190e27854706f770b5d7bfe68ea73"
PRE_OPTIMISATION_RATES_SHA256 = "fb6b7a124ba5d94f3456819f1987a2482017096b874d05206b0c9503e4413242"
OPTIMISED_RATES_SHA256 = "ccd3f274c02441c51729a069af645019f58d9c8a5457ba0a4a67e242e6c146b5"
# The BB170 boundary extension leaves every old positive-support draw intact.
# Its bridge applies only to requests that the earlier sampler could accept.
PRE_BOUNDARY_POPULATION_SHA256 = "1c3c24d76841a2bc770093c9c6ed412256f406c713407a5ef81eb9bd80b2fb71"
POSITIVE_LIFETIME_POPULATION_SHA256 = "cbb707ec8f375ed4382cb7402d09e9e89e3cdd2cb57e94efb0f818407df9422a"
THREAD_LIMITS = {name: "1" for name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"
)}
METRICS = ("bacteria", "sensed_mean", "expression_mean", "luminescence_mean",
           "luminescence_total", "free_signals", "bacteria_bound_signals",
           "polymer_bound_signals", "attached_polymers", "bridges", "largest_cluster")


def utc_now() -> str:
    """UTC labels provide chronology; elapsed time uses a monotonic clock."""
    return datetime.now(timezone.utc).isoformat()


@contextmanager
def controller_lock(directory: Path):
    """Allow only one controller per batch; POSIX releases the lock after a crash."""
    import fcntl

    with (directory / ".controller.lock").open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("This batch already has an active controller.") from error
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def available_cpus() -> int:
    """Respect process affinity and a SLURM per-task allocation when present."""
    count = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (
        os.cpu_count() or 1
    )
    allocation = os.environ.get("SLURM_CPUS_PER_TASK", "")
    if allocation.isdigit() and int(allocation) > 0:
        count = min(count, int(allocation))
    return max(1, count)


def validate_parameters(supplied: dict) -> dict:
    """Expand saved defaults and apply the same full-mode dataclass checks as CLI.

Batch files require the independent placement scheme explicitly through the
expanded defaults; they do not infer legacy placement from an old settings file.
"""
    classes = (InitialConfig, BindingConfig, MovementConfig, ResponseConfig,
               PopulationConfig, PlotConfig)
    settings = {"mode": "full"}
    for cls in classes:
        settings.update(asdict(cls()))
    unknown = set(supplied) - set(settings)
    if unknown:
        raise ValueError(f"Unknown model settings: {', '.join(sorted(unknown))}")
    settings.update(supplied)
    if settings["mode"] != "full" or settings["rng_scheme"] != "independent_v1":
        raise ValueError("Batches require mode=full and rng_scheme=independent_v1.")
    instances = [cls(**{k: settings[k] for k in asdict(cls())}) for cls in classes]
    initial, binding, _, _, population, plotting = instances
    binding.validate_polymer_rates(initial.polymers)
    if initial.bacteria > population.population_cap:
        raise ValueError("Initial bacteria exceed population_cap.")
    # Resolve fonts before spending time on any simulation. Store absolute font
    # paths so workers and subsequent resumes use the same supplied file.
    resolve_font(plotting)
    if settings["font_path"] is not None:
        settings["font_path"] = str(Path(settings["font_path"]).expanduser().resolve())
    return settings


def make_plan(spec: dict) -> list[dict]:
    """Create stable job identities independent of scheduling or worker count."""
    if not isinstance(spec, dict) or set(spec) != {"base_parameters", "conditions", "seeds"}:
        raise ValueError("Batch JSON requires exactly base_parameters, conditions, and seeds.")
    base, conditions, seeds = (spec[k] for k in ("base_parameters", "conditions", "seeds"))
    if not isinstance(base, dict) or "seed" in base:
        raise ValueError("base_parameters must be an object; put seeds in the seeds list.")
    if (not isinstance(seeds, list) or not seeds
            or any(type(s) is not int or not 0 <= s < 2**64 for s in seeds)
            or len(set(seeds)) != len(seeds)):
        raise ValueError("seeds must be a non-empty list of distinct integers in [0, 2**64).")
    if not isinstance(conditions, list) or not conditions:
        raise ValueError("conditions must be a non-empty list.")
    names, jobs = set(), []
    for condition in conditions:
        if not isinstance(condition, dict) or set(condition) != {"name", "overrides"}:
            raise ValueError("Each condition requires name and overrides.")
        name, overrides = condition["name"], condition["overrides"]
        if (not isinstance(name, str) or not re.fullmatch(r"[a-z][a-z0-9_]{0,39}", name)
                or name in names):
            raise ValueError("Condition names must be unique lowercase names (a-z, 0-9, _).")
        if not isinstance(overrides, dict) or "seed" in overrides:
            raise ValueError("Condition overrides must be an object without seed.")
        names.add(name)
        for seed in seeds:
            parameters = validate_parameters({**base, **overrides, "seed": seed})
            jobs.append({"id": f"{name}/seed_{seed}", "condition": name, "seed": seed,
                         "parameters": parameters, "status": "pending", "attempts": []})
    # Pointwise means and paired differences must use the same observation grid.
    grids = {(j["parameters"]["duration"], j["parameters"]["sample_interval"]) for j in jobs}
    if len(grids) != 1:
        raise ValueError("All conditions must use the same duration and sample_interval.")
    return jobs


def plan_hash(jobs: list[dict]) -> str:
    """Identify immutable requests separately from mutable execution records."""
    payload = [{k: j[k] for k in ("id", "condition", "seed", "parameters")} for j in jobs]
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def read_table(path: Path) -> list[dict[str, float]]:
    """Reject malformed or non-finite numerical observations rather than averaging them."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(handle)]
    if not rows or any(not math.isfinite(v) or v < 0 for row in rows for v in row.values()):
        raise ValueError(f"Empty, negative, or non-finite trajectory: {path}")
    return rows


def comparable_parameters(parameters: dict) -> dict:
    """Old runs implicitly used cutoff 10; only completed-run event budgets may differ."""
    return {"ps_crowding_threshold": CROWDING_THRESHOLD,
            **{k: v for k, v in parameters.items() if k != "max_events"}}


def numerical_sources_match(previous: dict, current: dict, parameters: dict) -> bool:
    """Accept exact sources or reviewed numerical-preserving predecessors only."""
    old, new = previous["source_sha256"], current["source_sha256"]
    if previous["numpy"] != current["numpy"]:
        return False
    for name in NUMERICAL_MODULES:
        if old.get(name) == new[name]:
            continue
        if name == "binding.py" and new[name] == PARAMETERISED_BINDING_SHA256:
            if old.get(name) == PRE_PROGRESS_BINDING_SHA256:
                continue
            if (old.get(name) == LEGACY_BINDING_SHA256
                    and parameters.get("ps_crowding_threshold", CROWDING_THRESHOLD) == CROWDING_THRESHOLD):
                continue
        if (name == "movement.py"
                and old.get(name) in (PRE_PROGRESS_MOVEMENT_SHA256, PROGRESS_MOVEMENT_SHA256)
                and new[name] == OPTIMISED_MOVEMENT_SHA256):
            continue
        if (name == "rates.py" and old.get(name) == PRE_OPTIMISATION_RATES_SHA256
                and new[name] == OPTIMISED_RATES_SHA256):
            continue
        if (name == "population.py" and old.get(name) == PRE_BOUNDARY_POPULATION_SHA256
                and new[name] == POSITIVE_LIFETIME_POPULATION_SHA256
                and parameters.get("doubling_time", 15000.0)
                > parameters.get("division_width", 20000.0) / 2):
            continue
        return False
    return True


def inspect_run(path: Path, parameters: dict, software: dict) -> tuple[dict, list[dict]]:
    """Check provenance, completion, observation grids, and signal/polymer balances.

    This is a boundary check on recorded results; the numerical solver performs
    the more detailed geometry and receptor-capacity checks during simulation.
    """
    meta = json.loads((path / "run.json").read_text(encoding="utf-8"))
    # A higher event budget cannot change a trajectory that already completed.
    # Keep the original budget in the copied run metadata for full provenance.
    if (meta.get("status") != "complete"
            or comparable_parameters(meta.get("parameters", {})) != comparable_parameters(parameters)):
        raise ValueError(f"Run is incomplete or its parameters differ: {path}")
    previous = meta["software"]
    if not numerical_sources_match(previous, software, parameters):
        raise ValueError(f"Numerical source or NumPy version differs: {path}")
    duration, interval = parameters["duration"], parameters["sample_interval"]
    times = [i * interval for i in range(math.floor(duration / interval) + 1)]
    if times[-1] != duration:
        times.append(duration)
    response, binding, population = [read_table(path / n) for n in (
        "response_trajectory.csv", "binding_trajectory.csv", "population_trajectory.csv"
    )]
    if any([r["time"] for r in table] != times for table in (response, binding, population)):
        raise ValueError(f"Missing, duplicated, or inconsistent sample times: {path}")
    initial_signals = 10000 if parameters["strain"] == "MM32" else 0
    merged = []
    for r, b, p in zip(response, binding, population):
        if not (
            b["free_signals"] + b["bacteria_bound_signals"] + b["polymer_bound_signals"]
            == b["total_signals"] == p["actual_total_signals"]
            == p["expected_total_signals"] == initial_signals + p["produced_signals"]
            and b["free_polymers"] + b["attached_polymers"] == parameters["polymers"]
            and b["pb_bonds"] == b["attached_polymers"] + b["bridges"]
            and r["sensed_total"] == b["bacteria_bound_signals"]
            and r["bacteria"] == p["bacteria"] == parameters["bacteria"] + p["divisions"]
            and math.isclose(r["luminescence_total"], r["luminescence_mean"] * r["bacteria"],
                             rel_tol=1e-12, abs_tol=1e-10)
            and math.isclose(r["sensed_total"], r["sensed_mean"] * r["bacteria"],
                             rel_tol=1e-12, abs_tol=1e-10)
            and math.isclose(r["expression_total"], r["expression_mean"] * r["bacteria"],
                             rel_tol=1e-12, abs_tol=1e-10)
        ):
            raise ValueError(f"Count/response balance failed at t={r['time']}: {path}")
        merged.append({**b, **p, **r})
    if meta["binding"]["stop_time"] != duration or meta["response"] != response[-1]:
        raise ValueError(f"Final metadata and observations disagree: {path}")
    return meta, merged


def reuse_runs(paths: list[Path], jobs: list[dict], software: dict) -> list[tuple[Path, dict]]:
    """Resolve every requested import before creating or launching a batch."""
    resolved, used = [], set()
    for path in paths:
        path = path.expanduser().resolve()
        parameters = json.loads((path / "run.json").read_text(encoding="utf-8"))["parameters"]
        matches = [j for j in jobs if comparable_parameters(j["parameters"])
                   == comparable_parameters(parameters)]
        if len(matches) != 1 or matches[0]["id"] in used:
            raise ValueError(f"Reused run must match exactly one unused planned job: {path}")
        inspect_run(path, matches[0]["parameters"], software)
        used.add(matches[0]["id"])
        resolved.append((path, matches[0]))
    return resolved


def matching_batch_runs(paths: list[Path], jobs: list[dict]) -> list[Path]:
    """Select compatible complete runs from previous batches by their settings.

    A previous condition's name is only a label. Match numerical/presentation
    requests instead, then let reuse_runs verify each original run's contents.
    The event budget may differ for a trajectory that already completed.
    """
    selected = []
    for directory in paths:
        directory = directory.expanduser().resolve()
        previous = json.loads((directory / "batch.json").read_text(encoding="utf-8"))
        if (previous.get("schema_version") != 1 or previous.get("status") != "complete"
                or previous.get("plan_sha256") != plan_hash(previous["jobs"])):
            raise ValueError(f"Reuse requires a complete batch with an unchanged plan: {directory}")
        matches = []
        for old in previous["jobs"]:
            parameters = old["parameters"]
            if old["status"] != "complete" or not any(
                comparable_parameters(j["parameters"]) == comparable_parameters(parameters)
                for j in jobs
            ):
                continue
            run = (directory / old["run_directory"]).resolve()
            if not run.is_relative_to(directory):
                raise ValueError(f"A reused run points outside its batch: {directory}")
            matches.append(run)
        if not matches:
            raise ValueError(f"No complete runs match this new plan: {directory}")
        print(f"Previous batch: {len(matches)} matching runs in {directory}")
        selected.extend(matches)
    return selected


def execute_jobs(directory: Path, manifest: dict, workers: int, progress="auto") -> int:
    """Bound concurrent child processes, record each attempt, and stop children on Ctrl-C."""
    pending = iter(j for j in manifest["jobs"] if j["status"] != "complete")
    active = []
    environment = {**os.environ, **THREAD_LIMITS}
    # Make source-tree use work from any launch directory, including an imported
    # API caller. Installed packages remain importable through the same Python.
    package_root = str(Path(__file__).resolve().parent.parent)
    environment["PYTHONPATH"] = package_root + os.pathsep + environment.get("PYTHONPATH", "")
    exhausted = False
    started = time.perf_counter()
    display = BatchProgress(directory, manifest["jobs"], progress)
    display.update(force=True)
    try:
        while active or not exhausted:
            while len(active) < workers and not exhausted:
                job = next(pending, None)
                if job is None:
                    exhausted = True
                    break
                number = len(job["attempts"]) + 1
                attempt = directory / job["id"] / f"attempt_{number:03d}"
                # A hard termination can leave a directory before its manifest
                # record is flushed. Preserve that attempt instead of overwriting it.
                while attempt.exists():
                    number += 1
                    attempt = directory / job["id"] / f"attempt_{number:03d}"
                attempt.mkdir(parents=True)
                write_json(attempt / "parameters.json", job["parameters"])
                record = {"directory": str(attempt.relative_to(directory)), "started_utc": utc_now()}
                job["attempts"].append(record)
                job["status"] = "running"
                log = (attempt / "console.log").open("w", encoding="utf-8")
                command = [sys.executable, "-m", "polymer_clustering", "--config",
                           str(attempt / "parameters.json"), "--output", str(attempt / "runs"),
                           "--progress", "off"]
                if display.console.enabled:
                    command += ["--progress-file", str(attempt / "progress.json")]
                try:
                    process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                               env=environment)
                except BaseException:
                    log.close()
                    job["status"] = "failed"
                    raise
                active.append((process, log, job, attempt, record))
                write_json(directory / "batch.json", manifest)
                display.message(f"Started {job['id']}")
            for item in list(active):
                process, log, job, attempt, record = item
                code = process.poll()
                if code is None:
                    continue
                log.close()
                record.update(exit_code=code, finished_utc=utc_now())
                runs = sorted((attempt / "runs").glob("*/run.json"))
                job["status"] = "failed"
                if len(runs) == 1:
                    run = runs[0].parent
                    job["run_directory"] = str(run.relative_to(directory))
                    try:
                        meta = json.loads(runs[0].read_text())
                        if code == 0:
                            inspect_run(run, job["parameters"], manifest["software"])
                            job["status"] = "complete"
                        elif code == 3 and meta["status"] == "incomplete":
                            job["status"] = "incomplete"
                        record["run_status"] = meta["status"]
                    except (OSError, ValueError, KeyError, TypeError) as error:
                        record["error"] = str(error)
                active.remove(item)
                write_json(directory / "batch.json", manifest)
                display.message(f"{job['status'].upper()}: {job['id']}; "
                                f"log: {attempt / 'console.log'}")
            display.update()
            if active:
                time.sleep(0.2)
    finally:
        # Do not leave expensive simulations running after the controller exits.
        for process, _, _, _, _ in active:
            if process.poll() is None:
                process.terminate()
        for process, log, job, _, record in active:
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            log.close()
            job["status"] = "interrupted"
            record.update(exit_code=process.returncode, finished_utc=utc_now())
        manifest["last_execution_wall_seconds"] = time.perf_counter() - started
        write_json(directory / "batch.json", manifest)
        display.close()
    statuses = {j["status"] for j in manifest["jobs"]}
    return 0 if statuses == {"complete"} else 3 if statuses <= {"complete", "incomplete"} else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--config", type=Path, help="Batch JSON: base_parameters, conditions, seeds.")
    source.add_argument("--resume", type=Path, help="Existing batch folder; replay unfinished jobs.")
    parser.add_argument("--workers", type=int, default=1, help="Concurrent runs; default 1.")
    parser.add_argument("--output", type=Path, default=Path("results/exercise5/replicates"))
    parser.add_argument("--reuse-run", type=Path, action="append", default=[],
                        help="Import a matching complete run folder; repeat for several runs.")
    parser.add_argument("--reuse-batch", type=Path, action="append", default=[],
                        help="Import matching runs from a complete batch folder; may be repeated.")
    parser.add_argument("--dry-run", action="store_true", help="Validate and list jobs; create no files.")
    add_progress_argument(parser)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Plan, launch/reuse, verify, and summarise an explicitly defined ensemble."""
    parser = build_parser()
    args = parser.parse_args(argv)
    overall_start = time.perf_counter()
    try:
        if not 1 <= args.workers <= available_cpus():
            raise ValueError(f"workers must be between 1 and the available CPU limit {available_cpus()}.")
        software = software_metadata()
        if args.resume:
            if args.reuse_run or args.reuse_batch:
                raise ValueError("Reuse flags apply to a new batch, not --resume.")
            directory = args.resume.expanduser().resolve()
            manifest = json.loads((directory / "batch.json").read_text())
            jobs = manifest["jobs"]
            if (manifest["plan_sha256"] != plan_hash(jobs)
                    or manifest["software"] != software):
                raise ValueError("Batch plan or software changed; use the original package to resume.")
            for job in jobs:
                validate_parameters(job["parameters"])
                if job["status"] == "complete":
                    inspect_run(directory / job["run_directory"], job["parameters"], software)
            imports = []
        else:
            spec = json.loads(args.config.read_text(encoding="utf-8"))
            jobs = make_plan(spec)
            prior_runs = matching_batch_runs(args.reuse_batch, jobs)
            imports = reuse_runs([*args.reuse_run, *prior_runs], jobs, software)
            if any(args.output.expanduser().resolve().is_relative_to(source)
                   for source, _ in imports):
                raise ValueError("The batch output parent cannot be inside an imported run.")
            manifest = {"schema_version": 1, "status": "planned", "started_utc": utc_now(),
                        "software": software, "plan_sha256": plan_hash(jobs), "jobs": jobs,
                        "child_thread_limits": THREAD_LIMITS,
                        "statistics": "Pointwise equal-weight means; sample SD (ddof=1); SEM=SD/sqrt(n). "
                                      "SD describes across-seed variation, not a confidence interval.",
                        "pairing": "Same seeds across conditions; seeds are independent within each "
                                   "condition. First condition is the paired-difference reference.",
                        "reuse_source_batches": [str(p.expanduser().resolve()) for p in args.reuse_batch],
                        "invocations": []}
        print(f"{len(jobs)} runs; workers={args.workers}; CPU limit={available_cpus()}.")
        imported_ids = {j["id"] for _, j in imports}
        for job in jobs:
            print(f"  {job['id']}: {'reuse' if job['id'] in imported_ids else job['status']}")
        retained = sum(j["id"] in imported_ids or j["status"] == "complete" for j in jobs)
        print(f"Ready: {retained} reused/retained; {len(jobs) - retained} simulations pending.")
        if args.dry_run:
            print("Plan validated; no simulations started and no output files created.")
            return 0
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.error(str(error))

    # All requests and imports are valid before allocating an output folder.
    try:
        if not args.resume:
            from .output import create_run_directory
            directory = create_run_directory(args.output, "full")
        with controller_lock(directory):
            if args.resume:
                # Refresh after acquiring the lock in case another controller
                # finished between our initial validation and lock acquisition.
                manifest = json.loads((directory / "batch.json").read_text())
                if manifest["plan_sha256"] != plan_hash(manifest["jobs"]) or (
                    manifest["software"] != software
                ):
                    raise ValueError("Batch plan or software changed before resume.")
                for job in manifest["jobs"]:
                    if job["status"] == "complete":
                        inspect_run(directory / job["run_directory"], job["parameters"], software)
            return finish_batch(directory, manifest, imports, args.workers, overall_start,
                                args.progress)
    except (OSError, ValueError, ImportError) as error:
        print(f"Could not start the batch: {error}", file=sys.stderr)
        return 1


def finish_batch(directory: Path, manifest: dict, imports: list, workers: int,
                 overall_start: float, progress="auto") -> int:
    """Persist the request before importing or launching; retain failures for resume."""
    from .batch_output import save_ensemble

    manifest["status"] = "running"
    manifest["invocations"].append({"started_utc": utc_now(), "workers": workers,
                                    "progress": progress})
    write_json(directory / "batch.json", manifest)
    print(f"Batch folder: {directory}", flush=True)
    code = 1
    try:
        for source, job in imports:
            destination = directory / job["id"] / "imported_run"
            shutil.copytree(source, destination)
            job.update(status="complete", run_directory=str(destination.relative_to(directory)),
                       reused_from=str(source))
            write_json(directory / "batch.json", manifest)
        code = execute_jobs(directory, manifest, workers, progress)
        if code == 0:
            print("All planned runs verified; saving ensemble tables and plots.", flush=True)
            save_ensemble(directory, manifest)
            manifest["status"] = "complete"
        else:
            manifest["status"] = "incomplete" if code == 3 else "failed"
            print("No ensemble average written: every planned run must complete.")
    except KeyboardInterrupt:
        code, manifest["status"] = 130, "interrupted"
        print("Batch interrupted; completed runs retained. Resume using the printed batch folder.")
    except Exception as error:
        manifest.update(status="failed", error=f"{type(error).__name__}: {error}")
        print(f"Batch failed: {error}", file=sys.stderr)
        code = 1
    finally:
        elapsed = time.perf_counter() - overall_start
        manifest["invocations"][-1].update(finished_utc=utc_now(), wall_seconds=elapsed)
        manifest["controller_wall_seconds_all_invocations"] = sum(
            i["wall_seconds"] for i in manifest["invocations"]
        )
        write_json(directory / "batch.json", manifest)
        print(f"Batch runtime this invocation: {elapsed:.3f} s. Outputs: {directory}")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
