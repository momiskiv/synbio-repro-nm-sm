# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Summarise independent runs without pooling cells or changing their light transform."""

import csv
from pathlib import Path

import matplotlib
import numpy as np

from .output import write_json
from .plotting import PlotConfig, plot_style, resolve_font


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write a complete table atomically; null variability becomes a blank CSV field."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def statistics(values: np.ndarray) -> dict:
    """Use sample SD across runs, leaving SD/SEM undefined for a single run."""
    n = len(values)
    sd = float(np.std(values, ddof=1)) if n > 1 else None
    return {"n": n, "mean": float(np.mean(values)), "sd": sd,
            "sem": sd / np.sqrt(n) if sd is not None else None}


def save_ensemble(directory: Path, manifest: dict) -> None:
    """Require the entire planned ensemble, then report per-condition and paired statistics."""
    # Import locally to keep the controller/renderer modules independently importable.
    from .batch import METRICS, inspect_run

    if any(j["status"] != "complete" for j in manifest["jobs"]):
        raise ValueError("Cannot average a batch with missing or incomplete runs.")
    grouped, records, final_rows = {}, {}, []
    for job in manifest["jobs"]:
        meta, rows = inspect_run(directory / job["run_directory"], job["parameters"],
                                 manifest["software"])
        grouped.setdefault(job["condition"], []).append(job)
        records[job["id"]] = rows
        peak = max(rows, key=lambda r: r["luminescence_mean"])
        final_rows.append({"condition": job["condition"], "seed": job["seed"],
                           "run_directory": job["run_directory"],
                           "reused": "reused_from" in job, "status": "complete",
                           "time": rows[-1]["time"],
                           **{k: rows[-1][k] for k in METRICS},
                           "sampled_peak_mean_light": peak["luminescence_mean"],
                           "sampled_peak_time": peak["time"],
                           "events": meta["binding"]["events"],
                           "run_wall_seconds": meta["timing"]["total_wall_seconds"],
                           "run_cpu_seconds": meta["timing"]["cpu_seconds"]})
    first = manifest["jobs"][0]
    times = [r["time"] for r in records[first["id"]]]
    if any([r["time"] for r in rows] != times for rows in records.values()):
        raise ValueError("Runs have different observation grids; no interpolation is performed.")
    ensemble, paired, pairwise = [], [], []
    reference = first["condition"]
    reference_jobs = {j["seed"]: j for j in grouped[reference]}
    for condition, jobs in grouped.items():
        for index, time in enumerate(times):
            for metric in METRICS:
                values = np.array([records[j["id"]][index][metric] for j in jobs])
                ensemble.append({"condition": condition, "time": time, "metric": metric,
                                 **statistics(values)})
                if condition != reference:
                    # Difference within each shared seed first; then summarise.
                    # Conditions sharing random streams are paired observations.
                    differences = np.array([
                        records[j["id"]][index][metric]
                        - records[reference_jobs[j["seed"]]["id"]][index][metric]
                        for j in jobs
                    ])
                    paired.append({"condition": condition, "reference": reference,
                                   "time": time, "metric": metric, **statistics(differences)})
    # With several affinity conditions, the mechanistic questions include P3
    # versus each disabled-binding control, not only versus the no-polymer case.
    # Emit each unordered pair once, always later-listed minus earlier-listed.
    # These are descriptive paired contrasts; no multiple-testing claim is made.
    names = list(grouped)
    if len(names) > 2:
        for number, condition in enumerate(names):
            for comparator in names[:number]:
                matched = {j["seed"]: j for j in grouped[comparator]}
                for index, time in enumerate(times):
                    for metric in METRICS:
                        differences = np.array([
                            records[j["id"]][index][metric]
                            - records[matched[j["seed"]]["id"]][index][metric]
                            for j in grouped[condition]
                        ])
                        pairwise.append({"condition": condition, "reference": comparator,
                                         "time": time, "metric": metric,
                                         **statistics(differences)})
    write_csv(directory / "replicate_summary.csv", final_rows)
    write_csv(directory / "ensemble_trajectory.csv", ensemble)
    if paired:
        write_csv(directory / "paired_differences.csv", paired)
    if pairwise:
        write_csv(directory / "pairwise_differences.csv", pairwise)
    plotting = {k: first["parameters"][k] for k in ("font_family", "font_path", "strict_font")}
    font = resolve_font(PlotConfig(**plotting))
    if font.fallback_used:
        print(font.notice())
    manifest["ensemble_plotting"] = font.metadata()
    # Original run CPU seconds include reused work. Keep that accounting
    # separate from the controller's actual wall time for this batch invocation.
    manifest["run_cpu_seconds_including_reused"] = sum(r["run_cpu_seconds"] for r in final_rows)
    manifest["newly_executed_run_cpu_seconds"] = sum(
        r["run_cpu_seconds"] for r in final_rows if not r["reused"]
    )
    manifest["aggregation"] = {"complete_runs": len(final_rows),
                               "replicates_per_condition": {k: len(v) for k, v in grouped.items()},
                               "paired_reference": reference,
                               "uncertainty": "Pointwise sample SD across independent seeds; "
                                              "not a confidence interval or across-cell SD.",
                               "weighting": "Equal weight per run; cell light transformed before "
                                            "within-run averaging; no baseline subtraction."}
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    colours = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
    affinity_style = {
        "no_polymer": ("No polymer", "#0072B2"),
        "p1_like_1000": ("P1-like (PB only)", "#009E73"),
        "p2_like_1000": ("P2-like (weak PS)", "#CC79A7"),
        "p3_no_pb_1000": ("P3 without PB", "#666666"),
        "p3_1000": ("P3: 1,000 polymers", "#D55E00"),
        "p3_2000": ("P3: 2,000 polymers", "#E69F00"),
        "p3_3000": ("P3: 3,000 polymers", "#009E73"),
        "p3_4000": ("P3: 4,000 polymers", "#56B4E9"),
        "p3_5000": ("P3: 5,000 polymers", "#8C6D31"),
        "p3_6000": ("P3: 6,000 polymers", "#666666"),
        "p3_7000": ("P3: 7,000 polymers", "#CC79A7"),
        "p3_cutoff_10": ("P3 (PS cutoff 10)", "#D55E00"),
        "p3_no_ps_cutoff": ("P3 (PS cutoff off)", "#0072B2"),
    }
    panels = (("luminescence_mean", "Mean light per cell", "Light (model a.u.)"),
              ("luminescence_total", "Total population light", "Light (model a.u.)"),
              ("sensed_mean", "Mean receptor occupancy", "Bound signals per cell"),
              ("largest_cluster", "Largest connected bacterial cluster", "Cells"))
    with plot_style(font):
        fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
        try:
            for ax, (metric, title, ylabel) in zip(axes.flat, panels):
                for number, (condition, jobs) in enumerate(grouped.items()):
                    selected = [r for r in ensemble
                                if r["condition"] == condition and r["metric"] == metric]
                    means = np.array([r["mean"] for r in selected])
                    label, colour = affinity_style.get(condition, (condition, colours[number % len(colours)]))
                    ax.plot(times, means, label=f"{label} (n={len(jobs)})", color=colour)
                    if len(jobs) > 1:
                        sd = np.array([r["sd"] for r in selected])
                        ax.fill_between(times, means - sd, means + sd, color=colour, alpha=0.15)
                ax.set(title=title, xlabel="Time (arbitrary units)", ylabel=ylabel,
                       xlim=(0, times[-1]))
                ax.grid(alpha=0.2)
                ax.legend(fontsize=8)
            changed_timing = any(j["parameters"]["doubling_time"] != 15000
                                 or j["parameters"]["division_width"] != 20000
                                 for j in manifest["jobs"])
            note = ("Division timing differs from source; diagnostic settings"
                    if changed_timing else "Raw model outputs; source reconstruction remains under review")
            if all(j["parameters"]["strain"] == "BB170"
                   and j["parameters"]["growth"]
                   and j["parameters"]["doubling_time"] == 10000
                   and j["parameters"]["division_width"] == 20000
                   for j in manifest["jobs"]):
                note = "BB170 faster growth: T=10000 a.u."
            if any(j["parameters"].get("ps_crowding_threshold", 10) != 10
                   for j in manifest["jobs"]):
                note = "PS crowding sensitivity: altered cutoff is an intervention, not a source default"
            uncertainty = ("means and pointwise ±1 SD" if any(len(jobs) > 1 for jobs in grouped.values())
                           else "means only (n=1)")
            fig.suptitle("Independent-run comparison | " + uncertainty + "\n" + note,
                         fontsize=13)
            fig.savefig(directory / "ensemble.png", dpi=200)
            fig.savefig(directory / "ensemble.svg")
        finally:
            plt.close(fig)
    write_json(directory / "batch.json", manifest)
