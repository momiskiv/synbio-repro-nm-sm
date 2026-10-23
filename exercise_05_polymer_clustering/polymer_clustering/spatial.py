# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 21 October 2026

"""Measure endpoint polymer crowding from completed runs, without simulating.

The original frame archives do not contain polymer coordinates over time.
Consequently these measurements describe the saved final state only. They can
identify spatial associations and occluded signal sites, but cannot establish
how much of a light difference was caused by crowding during the trajectory.
"""

import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time
from uuid import uuid4

import numpy as np

from . import __version__
from .batch import inspect_run, plan_hash
from .batch_output import statistics, write_csv
from .binding import BindingState, CROWDING_THRESHOLD, PS_PER_POLYMER, sample_counts
from .model import Bacterium
from .output import software_metadata, write_json
from .plotting import PlotConfig, resolve_font


STATE_FIELDS = ("occupancy", "free_signals", "bacteria_signals", "bacteria_polymers",
                "polymer_sites", "polymer_signals", "polymer_contacts", "bacteria_positions")
SOURCE_FILES = ("run.json", "binding_state.npz", "binding_trajectory.csv",
                "response_trajectory.csv", "population_trajectory.csv")
FINAL_METRICS = ("bacteria", "free_signals", "bacteria_bound_signals", "polymer_bound_signals",
                 "sensed_mean", "expression_mean", "luminescence_mean", "luminescence_total",
                 "attached_polymers", "bridges", "largest_cluster")


def file_hash(path: Path) -> str:
    """Hash inputs in bounded chunks, preserving provenance without loading large files."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fraction(numerator: float, denominator: float) -> float | None:
    """A zero denominator is undefined, not a measured zero fraction."""
    return numerator / denominator if denominator else None


def load_state(path: Path, parameters: dict, final: dict) -> BindingState:
    """Validate integer arrays, geometry, capacities and agreement with final counts."""
    with np.load(path, allow_pickle=False) as archive:
        if set(archive.files) != set(STATE_FIELDS):
            raise ValueError(f"Unexpected final-state array fields: {path}")
        arrays = {key: archive[key].copy() for key in STATE_FIELDS}
    if any(not np.issubdtype(a.dtype, np.integer) for a in arrays.values()):
        raise ValueError(f"Final-state arrays must contain integers: {path}")
    size, count = parameters["grid_size"], parameters["polymers"]
    if (arrays["occupancy"].shape != (size, size)
            or arrays["polymer_sites"].shape != (count,)
            or arrays["bacteria_positions"].shape != (int(final["bacteria"]), 3)):
        raise ValueError(f"Final-state dimensions disagree with recorded settings: {path}")
    positions = arrays.pop("bacteria_positions")
    identifiers = positions[:, 0].tolist()
    if len(set(identifiers)) != len(identifiers) or any(i < 0 for i in identifiers):
        raise ValueError(f"Invalid bacterial identities: {path}")
    state = BindingState(bacteria=tuple(Bacterium(*map(int, row)) for row in positions), **arrays)
    state.validate(int(final["total_signals"]), count)
    # Reuse the model's graph/count definitions; the archive and trajectory must
    # describe the same endpoint, including bridges and the largest component.
    observed = sample_counts(state, final["time"])
    if any(observed[key] != final[key] for key in observed):
        raise ValueError(f"Final arrays and final trajectory disagree: {path}")
    return state


def measure_state(state: BindingState, parameters: dict) -> tuple[dict, dict]:
    """Return endpoint metrics and auditable site arrays using existing model rules.

    PB and PS occupy separate capacities. Empty PS slots become inaccessible
    at the run's cutoff (default P >= 10); already-bound signals retain their off rate.
    The hypothetical propensity without this cutoff holds the entire state
    fixed. It is a local rate calculation, never a simulated intervention.
    """
    shape, sites = state.occupancy.shape, state.polymer_sites
    density = np.bincount(sites, minlength=state.occupancy.size).reshape(shape)
    threshold = parameters.get("ps_crowding_threshold", CROWDING_THRESHOLD)
    crowded = density >= threshold if threshold is not None else np.zeros(shape, dtype=bool)
    source_dense = density >= CROWDING_THRESHOLD
    on_bacteria = state.occupancy >= 0
    free_slots = PS_PER_POLYMER - state.polymer_signals
    empty = np.bincount(sites, weights=free_slots, minlength=density.size).reshape(shape)
    bound = np.bincount(sites, weights=state.polymer_signals,
                        minlength=density.size).reshape(shape)
    blocked = np.where(crowded, empty, 0)
    accessible = empty - blocked
    on_rate, off_rate = (parameters.get(key) or 0.0 for key in ("k_ps_on", "k_ps_off"))
    unblocked_rate = on_rate * empty * state.free_signals
    binding_rate = on_rate * accessible * state.free_signals
    npolymers = len(sites)
    crowded_polymers = int(density[crowded].sum())
    bacterial_polymers = int(density[on_bacteria].sum())
    values = {
        "polymers": npolymers,
        "max_local_polymers": int(density.max()),
        "polymer_occupied_sites": int(np.count_nonzero(density)),
        "crowded_sites": int(crowded.sum()),
        "polymers_at_crowded_sites": crowded_polymers,
        "polymer_crowded_fraction": fraction(crowded_polymers, npolymers),
        # A fixed reference density remains meaningful when the active cutoff
        # is disabled: disabling a rule does not make polymers less dense.
        "sites_at_source_density_threshold": int(source_dense.sum()),
        "polymer_dense_fraction_at_10": fraction(int(density[source_dense].sum()), npolymers),
        "polymers_on_bacterial_sites": bacterial_polymers,
        "polymer_on_bacteria_fraction": fraction(bacterial_polymers, npolymers),
        "empty_ps_slots": int(empty.sum()),
        "blocked_empty_ps_slots": int(blocked.sum()),
        "blocked_empty_ps_fraction": fraction(float(blocked.sum()), float(empty.sum())),
        "ps_capacity_occupancy": fraction(int(state.polymer_signals.sum()),
                                           PS_PER_POLYMER * npolymers),
        "ps_bound_on_bacterial_sites": int(bound[on_bacteria].sum()),
        "ps_bind_propensity": float(binding_rate.sum()),
        "ps_bind_propensity_without_cutoff": float(unblocked_rate.sum()),
        "ps_bind_propensity_lost_to_cutoff": float((unblocked_rate - binding_rate).sum()),
        "ps_unbind_propensity": float(off_rate * state.polymer_signals.sum()),
        "mean_free_signals_bacterial_sites": (
            float(state.free_signals[on_bacteria].mean()) if on_bacteria.any() else None),
        "mean_free_signals_other_sites": (
            float(state.free_signals[~on_bacteria].mean()) if (~on_bacteria).any() else None),
    }
    maps = {"occupancy": state.occupancy, "free_signals": state.free_signals,
            "polymer_density": density, "polymer_bound_signals": bound,
            "empty_ps_slots": empty, "blocked_empty_ps_slots": blocked,
            "ps_bind_propensity": binding_rate,
            "ps_bind_propensity_without_cutoff": unblocked_rate}
    return values, maps


def select_jobs(manifest: dict, conditions: list[str] | None,
                seeds: list[int] | None) -> list[dict]:
    """Require a complete, unchanged plan and balanced, uniquely identified seed pairs."""
    jobs = manifest.get("jobs", [])
    if (manifest.get("schema_version") != 1 or manifest.get("status") != "complete"
            or not jobs or any(j.get("status") != "complete" for j in jobs)
            or manifest.get("plan_sha256") != plan_hash(jobs)):
        raise ValueError("Spatial analysis requires a complete batch with an unchanged plan.")
    if (len({j["id"] for j in jobs}) != len(jobs)
            or len({(j["condition"], j["seed"]) for j in jobs}) != len(jobs)
            or any(j["seed"] != j["parameters"]["seed"] for j in jobs)):
        raise ValueError("Duplicate or inconsistent condition/seed records.")
    available = list(dict.fromkeys(j["condition"] for j in jobs))
    names = conditions if conditions is not None else available
    if not names or len(set(names)) != len(names) or set(names) - set(available):
        raise ValueError(f"Choose distinct conditions from: {', '.join(available)}")
    requested = set(seeds) if seeds is not None else None
    if seeds is not None and (not seeds or len(requested) != len(seeds)):
        raise ValueError("Choose a non-empty list of distinct seeds.")
    selected = [j for name in names for j in jobs if j["condition"] == name
                and (requested is None or j["seed"] in requested)]
    sets = [{j["seed"] for j in selected if j["condition"] == name} for name in names]
    if (not sets[0] or any(s != sets[0] for s in sets)
            or requested is not None and sets[0] != requested):
        raise ValueError("Selected conditions must contain the same requested seed set.")
    settings = {(j["parameters"]["duration"], j["parameters"]["grid_size"],
                 j["parameters"]["strain"]) for j in selected}
    if len(settings) != 1:
        raise ValueError("Endpoint comparisons require identical duration, grid size and strain.")
    return selected


def aggregate(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Weight runs equally and compute paired differences before their sample SD.

    Blank ratios stay undefined, including the no-polymer denominator. Their
    aggregate n is zero, while actual zero counts remain valid observations.
    """
    names = list(dict.fromkeys(row["condition"] for row in rows))
    metrics = [k for k in rows[0] if k not in ("condition", "seed", "run_directory", "time")]
    grouped = {name: {r["seed"]: r for r in rows if r["condition"] == name} for name in names}

    def describe(values: list[float]) -> dict:
        return statistics(np.array(values)) if values else {
            "n": 0, "mean": None, "sd": None, "sem": None}

    ensemble = [{"condition": name, "time": rows[0]["time"], "metric": metric,
                 **describe([r[metric] for r in grouped[name].values() if r[metric] is not None])}
                for name in names for metric in metrics]
    paired = []
    for earlier, later in combinations(names, 2):
        for metric in metrics:
            values = [grouped[later][seed][metric] - row[metric]
                      for seed, row in grouped[earlier].items()
                      if row[metric] is not None and grouped[later][seed][metric] is not None]
            paired.append({"condition": later, "reference": earlier, "time": rows[0]["time"],
                           "metric": metric, **describe(values)})
    return ensemble, paired


def analyse(batch: Path, conditions: list[str] | None, seeds: list[int] | None,
            map_seed: int | None) -> tuple[list[dict], list[dict], list[dict], dict]:
    """Read and validate every selected source before creating an output directory."""
    manifest = json.loads((batch / "batch.json").read_text(encoding="utf-8"))
    jobs = select_jobs(manifest, conditions, seeds)
    if map_seed is None:
        map_seed = min(j["seed"] for j in jobs)
    if map_seed not in {j["seed"] for j in jobs}:
        raise ValueError("The map seed must be one of the selected seeds.")
    software = software_metadata()
    rows, site_rows, maps, provenance = [], [], [], []
    for job in jobs:
        directory = (batch / job["run_directory"]).resolve()
        if not directory.is_relative_to(batch) or directory == batch:
            raise ValueError("A recorded run path escapes the batch directory.")
        # Resolve individual files too: a symlink must not silently substitute
        # data from outside the selected batch.
        if any(not (directory / name).resolve().is_relative_to(batch) for name in SOURCE_FILES):
            raise ValueError("A recorded input file escapes the batch directory.")
        metadata, history = inspect_run(directory, job["parameters"], software)
        final = history[-1]
        state = load_state(directory / "binding_state.npz", job["parameters"], final)
        values, arrays = measure_state(state, job["parameters"])
        rows.append({"condition": job["condition"], "seed": job["seed"],
                     "run_directory": job["run_directory"], "time": final["time"],
                     **{k: final[k] for k in FINAL_METRICS}, **values})
        for y, x in np.ndindex(state.occupancy.shape):
            site_rows.append({"condition": job["condition"], "seed": job["seed"],
                              "time": final["time"], "x": x, "y": y,
                              **{key: array[y, x].item() for key, array in arrays.items()}})
        if job["seed"] == map_seed:
            maps.append({"condition": job["condition"], "seed": map_seed,
                         "time": final["time"], "arrays": arrays,
                         "ps_crowding_threshold": job["parameters"].get(
                             "ps_crowding_threshold", CROWDING_THRESHOLD),
                         "bacteria": state.bacteria})
        provenance.append({"condition": job["condition"], "seed": job["seed"],
                           "run_directory": job["run_directory"],
                           "source_program_version": metadata["software"]["program_version"],
                           "parameters": job["parameters"],
                           "input_sha256": {name: file_hash(directory / name)
                                            for name in SOURCE_FILES}})
    report = {"schema_version": 1, "status": "validated", "analysis_version": __version__,
              "scope": "Final-state spatial diagnostics; no simulation executed.",
              "source_batch": str(batch), "source_batch_sha256": file_hash(batch / "batch.json"),
              "source_plan_sha256": manifest["plan_sha256"], "map_seed": map_seed,
              "software": software, "inputs": provenance,
              "rules": {"crowding_threshold": CROWDING_THRESHOLD,
                        "active_threshold": "Per run: inputs[].parameters.ps_crowding_threshold; "
                        "missing means 10, null disables PS occlusion.",
                        "ps_slots_per_polymer": PS_PER_POLYMER,
                        "on_bacteria": "Polymer anchor lies on a bacterial footprint.",
                        "blocked_empty_ps_slots": "Unoccupied PS slots blocked by the active "
                        "per-run cutoff; zero when disabled.",
                        "without_cutoff_propensity": "Hypothetical PS on-rate at this same fixed "
                        "state with all unoccupied PS slots accessible; not a rerun.",
                        "uncertainty": "Equal-weight run means and sample SD (ddof=1), not CI. "
                        "Paired differences are later-listed minus earlier-listed conditions.",
                        "missing": "Zero denominators yield blank ratios, with n=0 when all "
                        "values are undefined; one replicate has blank SD/SEM."},
              "limitations": ["No polymer time course is inferred from the final state.",
                              "Spatial associations and instantaneous propensities do not "
                              "isolate crowding's causal contribution to light suppression.",
                              "No physical concentration or experimental light calibration."]}
    return rows, site_rows, maps, report


def main(argv: list[str] | None = None) -> int:
    """Expose post-processing separately from simulation and retain input files intact."""
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--batch", type=Path, required=True, help="Completed folder containing batch.json.")
    parser.add_argument("--conditions", nargs="+", help="Condition names; default all, in batch order.")
    parser.add_argument("--seeds", nargs="+", type=int, help="Matched seeds; default all.")
    parser.add_argument("--map-seed", type=int, help="Illustrative map seed; default smallest selected seed.")
    parser.add_argument("--output", type=Path, default=Path("results/exercise5/spatial_diagnostics"))
    parser.add_argument("--font-family", default="Arial")
    parser.add_argument("--font-path")
    parser.add_argument("--strict-font", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs; write nothing.")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args(argv)
    started = time.perf_counter()
    cpu_started = time.process_time()
    directory = None
    try:
        batch, output = args.batch.expanduser().resolve(), args.output.expanduser().resolve()
        if output.is_relative_to(batch):
            raise ValueError("Choose an output root outside the original batch folder.")
        font = resolve_font(PlotConfig(args.font_family, args.font_path, args.strict_font))
        rows, site_rows, maps, report = analyse(batch, args.conditions, args.seeds, args.map_seed)
        if font.fallback_used:
            print(font.notice())
        print(f"Validated {len(rows)} completed runs; endpoint t={rows[0]['time']:g} a.u.")
        if args.dry_run:
            print("No simulations started and no output files created.")
            print(f"Validation time: {time.perf_counter() - started:.3f} s.")
            return 0
        ensemble, paired = aggregate(rows)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        directory = output / f"spatial_{stamp}_{uuid4().hex[:8]}"
        directory.mkdir(parents=True)
        write_csv(directory / "spatial_replicates.csv", rows)
        write_csv(directory / "spatial_ensemble.csv", ensemble)
        write_csv(directory / "spatial_sites.csv", site_rows)
        if paired:
            write_csv(directory / "spatial_paired_differences.csv", paired)
        from .spatial_output import save_plots, save_bundle
        save_plots(directory, rows, ensemble, maps, font)
        report.update(status="complete", plotting=font.metadata(),
                      timing={"wall_seconds": time.perf_counter() - started,
                              "cpu_seconds": time.process_time() - cpu_started,
                              "scope": "Validation, analysis and rendering; excludes final "
                              "report write, summary ZIP creation and terminal printing."})
        write_json(directory / "spatial_report.json", report)
        save_bundle(directory)
        print(f"Spatial results: {directory}")
        print(f"Upload: {directory / 'spatial_summary.zip'}")
        print(f"Runtime: {time.perf_counter() - started:.3f} s; no simulations executed.")
        return 0
    except (OSError, ValueError, KeyError, TypeError, IndexError, json.JSONDecodeError) as error:
        if directory is not None:
            write_json(directory / "spatial_report.json", {
                "status": "failed", "error": str(error), "analysis_version": __version__})
        print(f"Could not complete spatial analysis: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
