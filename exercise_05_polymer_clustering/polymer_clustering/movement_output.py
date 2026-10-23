# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Spatial tables, start/end maps, and a labelled animation of observed states."""

import csv
from pathlib import Path

import matplotlib

import numpy as np
from PIL import Image

from .plotting import styled_plot
from .movement import MovementResult
from .model import Bacterium
from .output import freeze_animation_layout


def save_movement(directory: Path, result: MovementResult) -> None:
    """Write movement outputs beside the common binding outputs in the run folder."""
    for filename, rows, fields in (
        ("movement_positions.csv", result.positions, ["time", "bacterium", "x", "y"]),
        ("final_bacteria.csv", [{"bacterium": b.identifier, "x": b.x, "y": b.y}
                                for b in result.state.bacteria], ["bacterium", "x", "y"]),
    ):
        with (directory / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    # Frames are observations, not time integration steps. Increasing their
    # number does not change the event sequence or physical elapsed time.
    np.savez_compressed(
        directory / "movement_frames.npz",
        times=np.array([frame.time for frame in result.frames]),
        free_signals=np.stack([frame.free_signals for frame in result.frames]),
        bacteria_positions=frame_positions(result),
        bacteria_counts=np.array([len(frame.bacteria) for frame in result.frames]),
    )
    save_spatial_views(directory, result)


def frame_positions(result: MovementResult) -> np.ndarray:
    """Pad unborn cells with -1, keeping a rectangular numeric frame archive."""
    positions = np.full((len(result.frames), len(result.state.bacteria), 3), -1, dtype=np.int64)
    for index, frame in enumerate(result.frames):
        positions[index, :len(frame.bacteria)] = np.array(
            [(b.identifier, b.x, b.y) for b in frame.bacteria], dtype=np.int64,
        ).reshape(-1, 3)
    return positions


@styled_plot
def save_spatial_views(directory: Path, result: MovementResult) -> None:
    """Draw comparable free-signal maps and animate the same sampled states."""
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.axes import Axes
    from matplotlib.patches import Rectangle

    size = result.state.occupancy.shape[0]
    maximum = max(1, max(int(frame.free_signals.max()) for frame in result.frames))
    suffix = " | INCOMPLETE" if result.status != "complete" else ""
    if result.scenario_note:
        suffix += "\n" + result.scenario_note

    def axes_setup(ax: Axes) -> None:
        ax.set(xlim=(0, size), ylim=(0, size), xlabel="x (lattice units)",
               ylabel="y (lattice units)")
        ax.set_aspect("equal")

    def cell_outline(ax: Axes, bacterium: Bacterium) -> Rectangle:
        # ID colours stay fixed between frames, making cell translations visible.
        color = plt.get_cmap("tab20")(bacterium.identifier % 20)
        rectangle = Rectangle((bacterium.x, bacterium.y), 2, 2,
                              facecolor="none", edgecolor=color, linewidth=1.5)
        ax.add_patch(rectangle)
        return rectangle

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), layout="constrained")
    try:
        for ax, frame, label in zip(axes, (result.frames[0], result.frames[-1]), ("Start", "End")):
            axes_setup(ax)
            image = ax.imshow(frame.free_signals, origin="lower", extent=(0, size, 0, size),
                              interpolation="nearest", cmap="cividis", vmin=0, vmax=maximum)
            for bacterium in frame.bacteria:
                cell_outline(ax, bacterium)
            ax.set_title(f"{label}: t = {frame.time:.3g} a.u.")
        fig.colorbar(image, ax=axes, label="Free signals per lattice site", shrink=0.8)
        fig.suptitle("Binding and movement | outlines: bacteria" + suffix)
        fig.savefig(directory / "movement.png", dpi=180)
        fig.savefig(directory / "movement.svg")
    finally:
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5.6), layout="constrained")
    pictures: list[Image.Image] = []
    try:
        axes_setup(ax)
        image = ax.imshow(result.frames[0].free_signals, origin="lower", extent=(0, size, 0, size),
                          interpolation="nearest", cmap="cividis", vmin=0, vmax=maximum)
        fig.colorbar(image, ax=ax, label="Free signals per lattice site", shrink=0.8)
        rectangles = [cell_outline(ax, bacterium) for bacterium in result.state.bacteria]
        for rectangle in rectangles:
            rectangle.set_visible(False)
        ax.set_title("Binding + movement | t = 0 a.u." + suffix)
        freeze_animation_layout(fig)
        for frame in result.frames:
            image.set_data(frame.free_signals)
            for index, rectangle in enumerate(rectangles):
                rectangle.set_visible(index < len(frame.bacteria))
            for rectangle, bacterium in zip(rectangles, frame.bacteria):
                rectangle.set_xy((bacterium.x, bacterium.y))
            ax.set_title(f"Binding + movement | t = {frame.time:.3g} a.u." + suffix)
            fig.canvas.draw()
            pixels = np.asarray(fig.canvas.buffer_rgba()).copy()
            pictures.append(Image.fromarray(pixels).convert("RGB"))
        # Playback speed is illustrative: the frame title records simulation
        # time. GIF duration is not a conversion from a.u. to experimental hours.
        pictures[0].save(directory / "movement.gif", save_all=True,
                         append_images=pictures[1:], duration=200, loop=0)
    finally:
        plt.close(fig)
        for picture in pictures:
            picture.close()
