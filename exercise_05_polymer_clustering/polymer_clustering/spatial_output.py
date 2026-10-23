# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 21 October 2026

"""Comparable endpoint maps, replicate plots and a compact review ZIP."""

from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import matplotlib
import numpy as np

from .binding import CROWDING_THRESHOLD
from .plotting import ResolvedFont, plot_style


STYLES = {
    "no_polymer": ("No polymer", "#0072B2"),
    "p1_like_1000": ("P1-like\nPB only", "#009E73"),
    "p2_like_1000": ("P2-like\nweak PS", "#CC79A7"),
    "p3_no_pb_1000": ("P3 without PB", "#666666"),
    "p3_1000": ("P3\nPB + PS", "#D55E00"),
    "p3_cutoff_10": ("P3\ncutoff 10", "#D55E00"),
    "p3_no_ps_cutoff": ("P3\nPS cutoff off", "#0072B2"),
}


def save_plots(directory: Path, rows: list[dict], ensemble: list[dict],
               maps: list[dict], font: ResolvedFont) -> None:
    """Use shared map colour limits and distinguish an example seed from the ensemble."""
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle

    names = list(dict.fromkeys(r["condition"] for r in rows))
    summary = {(r["condition"], r["metric"]): r for r in ensemble}
    labels = [STYLES.get(name, (name, "#0072B2"))[0] for name in names]
    colours = [STYLES.get(name, (name, "#0072B2"))[1] for name in names]
    panels = (("polymer_dense_fraction_at_10", "Polymers at P ≥ 10 per site", "Polymers (%)", 100),
              ("blocked_empty_ps_slots", "Empty PS slots blocked by crowding", "PS slots", 1),
              ("ps_capacity_occupancy", "Occupied polymer signal capacity", "PS capacity (%)", 100),
              ("luminescence_mean", "Mean light per cell", "Light (model a.u.)", 1))
    with plot_style(font):
        fig, axes = plt.subplots(2, 2, figsize=(max(9, len(names) * 2.1), 8), layout="constrained")
        try:
            for ax, (metric, title, ylabel, scale) in zip(axes.flat, panels):
                for index, (name, colour) in enumerate(zip(names, colours)):
                    stats = summary[name, metric]
                    if stats["mean"] is None:
                        ax.text(index, 0.02, "n/a", ha="center", va="bottom",
                                transform=ax.get_xaxis_transform())
                        continue
                    ax.bar(index, stats["mean"] * scale, color=colour, alpha=0.25, width=0.65)
                    if stats["sd"] is not None:
                        ax.errorbar(index, stats["mean"] * scale, yerr=stats["sd"] * scale,
                                    fmt="none", ecolor=colour, capsize=4)
                    values = [r[metric] * scale for r in rows
                              if r["condition"] == name and r[metric] is not None]
                    # Deterministic display offsets are not simulation random draws.
                    offsets = np.linspace(-0.13, 0.13, len(values))
                    ax.scatter(index + offsets, values, color=colour, s=20, zorder=3)
                ax.set(title=title, ylabel=ylabel, ylim=(0, None),
                       xlim=(-0.6, len(names) - 0.4))
                if all(summary[name, metric]["mean"] in (None, 0) for name in names):
                    # An all-zero endpoint should not acquire fractional count
                    # ticks merely from Matplotlib's default zero-range padding.
                    ax.set_ylim(0, 1)
                    ax.set_yticks([0, 1])
                ax.set_xticks(range(len(names)), labels)
                ax.grid(axis="y", alpha=0.2)
            fig.suptitle(f"Final-state comparison at t = {rows[0]['time']:g} a.u.\n"
                         "Each point: one seed; bars: means ± sample SD", fontsize=13)
            fig.savefig(directory / "spatial_comparison.png", dpi=180)
            fig.savefig(directory / "spatial_comparison.svg")
        finally:
            plt.close(fig)

        fields = (("polymer_density", "Polymer density", "Polymers per site", "viridis"),
                  ("blocked_empty_ps_slots", "Blocked empty PS slots", "Slots per site", "magma"),
                  ("free_signals", "Free signals", "Signals per site", "cividis"))
        # A column has one normalisation for every condition, preventing a dense
        # map from appearing deceptively similar to a sparse map through autoscaling.
        maxima = [max(1, max(float(m["arrays"][key].max()) for m in maps))
                  for key, _, _, _ in fields]
        fig, axes = plt.subplots(len(maps), 3, squeeze=False,
                                 figsize=(12, max(4.2, 3.6 * len(maps))), layout="constrained")
        try:
            for row, example in enumerate(maps):
                arrays = example["arrays"]
                size = arrays["occupancy"].shape[0]
                for column, (key, title, _, cmap) in enumerate(fields):
                    ax = axes[row, column]
                    images = ax.imshow(arrays[key], origin="lower", extent=(0, size, 0, size),
                                       interpolation="nearest", cmap=cmap, vmin=0, vmax=maxima[column])
                    for b in example["bacteria"]:
                        ax.add_patch(Rectangle((b.x, b.y), 2, 2, fill=False,
                                               edgecolor="white", linewidth=0.55))
                    ax.set(xlim=(0, size), ylim=(0, size), xlabel="x (lattice units)",
                           ylabel="y (lattice units)", aspect="equal")
                    label = STYLES.get(example["condition"], (example["condition"], ""))[0]
                    cutoff = example.get("ps_crowding_threshold", CROWDING_THRESHOLD)
                    cutoff_label = "off" if cutoff is None else str(cutoff)
                    ax.set_title(label.replace("\n", " ") + " | " + title, fontsize=10)
                    if column == 0:
                        ax.text(0.02, 0.02, f"PS cutoff: {cutoff_label}", transform=ax.transAxes,
                                color="white", fontsize=8,
                                bbox={"facecolor": "black", "alpha": .6, "edgecolor": "none"})
                    if row == len(maps) - 1:
                        fig.colorbar(images, ax=list(axes[:, column]), label=fields[column][2],
                                     shrink=0.8)
            fig.suptitle(f"Saved endpoint | seed {maps[0]['seed']} | t = {maps[0]['time']:g} a.u.\n"
                         "White outlines: bacteria. Active PS cutoff labelled on each row.",
                         fontsize=12)
            fig.savefig(directory / "spatial_maps.png", dpi=180)
            fig.savefig(directory / "spatial_maps.svg")
        finally:
            plt.close(fig)


def save_bundle(directory: Path) -> None:
    """Package compact evidence; keep the full site table and vector maps locally."""
    names = ("spatial_report.json", "spatial_replicates.csv", "spatial_ensemble.csv",
             "spatial_paired_differences.csv", "spatial_comparison.png", "spatial_maps.png")
    with ZipFile(directory / "spatial_summary.zip", "w", compression=ZIP_DEFLATED) as archive:
        for name in names:
            path = directory / name
            if path.is_file():
                archive.write(path, arcname=name)
