# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 21 October 2026

"""Save inspectable binding histories and final state in the current run folder."""

import csv
from pathlib import Path

import matplotlib

import numpy as np

from .plotting import styled_plot
from .binding import BindingResult


def save_binding(directory: Path, result: BindingResult, moving: bool = False,
                 scenario_note: str = "") -> None:
    """Persist sampled counts, individual events, polymer contacts, and arrays."""
    # Keep the pre-binding initial_state.npz intact; these are the final arrays.
    state = result.state
    np.savez_compressed(
        directory / "binding_state.npz", occupancy=state.occupancy,
        free_signals=state.free_signals, bacteria_signals=state.bacteria_signals,
        bacteria_polymers=state.bacteria_polymers, polymer_sites=state.polymer_sites,
        polymer_signals=state.polymer_signals, polymer_contacts=state.polymer_contacts,
        bacteria_positions=np.array([(b.identifier, b.x, b.y) for b in state.bacteria],
                                    dtype=np.int64).reshape(-1, 3),
    )
    event_fields = ["event", "time", "reaction", "x", "y", "polymer", "bacterium"]
    if moving:
        event_fields += ["to_x", "to_y", "cluster_size", "cluster_members",
                         "advected_signals", "advected_polymers", "daughter"]
    for filename, rows, fields in (
        ("binding_trajectory.csv", result.trajectory, list(result.trajectory[0])),
        ("binding_events.csv", result.events,
         event_fields),
    ):
        with (directory / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    # A row per polymer makes bridge identity reviewable without reading NPZ.
    size = state.occupancy.shape[0]
    with (directory / "polymer_binding.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["polymer", "x", "y", "bound_signals", "contact_0", "contact_1",
                  "bacterium_0", "bacterium_1"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for polymer, site in enumerate(state.polymer_sites):
            y, x = divmod(int(site), size)
            contacts = state.polymer_contacts[polymer]
            bacteria = [int(state.occupancy.ravel()[c]) if c >= 0 else -1 for c in contacts]
            writer.writerow(dict(zip(fields, [polymer, x, y, int(state.polymer_signals[polymer]),
                                               *contacts.tolist(), *bacteria])))
    with (directory / "bacteria_binding.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["bacterium", "bound_signals", "pb_bonds"])
        writer.writeheader()
        for bacterium in state.bacteria:
            footprint = np.s_[bacterium.y:bacterium.y + 2, bacterium.x:bacterium.x + 2]
            writer.writerow({"bacterium": bacterium.identifier,
                             "bound_signals": int(state.bacteria_signals[footprint].sum()),
                             "pb_bonds": int(state.bacteria_polymers[footprint].sum())})
    save_binding_plot(directory, result, moving, scenario_note)


@styled_plot
def save_binding_plot(directory: Path, result: BindingResult, moving: bool = False,
                      scenario_note: str = "") -> None:
    """Plot fixed-position diagnostic counts without implying luminescence."""
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    times = [row["time"] for row in result.trajectory]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    try:
        for key, label, color in (
            ("free_signals", "Free", "#0072B2"),
            ("bacteria_bound_signals", "Bacterium-bound", "#D55E00"),
            ("polymer_bound_signals", "Polymer-bound", "#009E73"),
        ):
            axes[0].step(times, [row[key] for row in result.trajectory], where="post",
                         label=label, color=color)
        for key, label, color in (
            ("attached_polymers", "Attached polymers", "#0072B2"),
            ("pb_bonds", "Bacterium-polymer bonds", "#D55E00"),
            ("bridges", "Bridging polymers", "#009E73"),
        ):
            axes[1].step(times, [row[key] for row in result.trajectory], where="post",
                         label=label, color=color)
        axes[0].set(title="Signal pools", ylabel="Signal count")
        axes[1].set(title="Polymer attachments", ylabel="Count")
        for ax in axes:
            ax.set_xlabel("Time (arbitrary units)")
            ax.set_ylim(bottom=0)
            ax.legend()
            ax.grid(alpha=0.2)
        suffix = " | INCOMPLETE: event limit" if result.status != "complete" else ""
        if scenario_note:
            suffix += "\n" + scenario_note
        scope = ("Binding and movement | signal pools and attachments" if moving
                 else "Binding only: fixed positions; no growth or light calculation")
        fig.suptitle(scope + suffix)
        fig.savefig(directory / "binding.png", dpi=180)
        fig.savefig(directory / "binding.svg")
    finally:
        plt.close(fig)
