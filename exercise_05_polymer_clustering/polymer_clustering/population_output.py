# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Population history, lineage, and an explicit signal-production balance."""

import csv
from pathlib import Path

import matplotlib
from .plotting import styled_plot
from .movement import MovementResult


@styled_plot
def save_population(directory: Path, result: MovementResult) -> None:
    """Make births and the changing signal conservation target easy to audit."""
    for name, rows, fields in (
        ("population_trajectory.csv", result.population_trajectory,
         list(result.population_trajectory[0])),
        ("lineage.csv", result.lineage, ["bacterium", "parent", "birth_time", "generation"]),
    ):
        with (directory / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    rows = result.population_trajectory
    times = [row["time"] for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), layout="constrained")
    try:
        axes[0].step(times, [row["bacteria"] for row in rows], where="post", color="#0072B2")
        axes[0].set(ylabel="Number of bacteria", title="Age-dependent division")
        axes[1].step(times, [row["produced_signals"] for row in rows], where="post", color="#009E73")
        axes[1].set(ylabel="Cumulative new signals", title="Signal production")
        for ax in axes:
            ax.set_xlabel("Time (arbitrary units)")
            ax.set_ylim(bottom=0)
            ax.grid(alpha=0.2)
        suffix = " | INCOMPLETE" if result.status != "complete" else ""
        if result.scenario_note:
            suffix += "\n" + result.scenario_note
        fig.suptitle("Population history | sampled observations" + suffix)
        fig.savefig(directory / "population.png", dpi=180)
        fig.savefig(directory / "population.svg")
    finally:
        plt.close(fig)
