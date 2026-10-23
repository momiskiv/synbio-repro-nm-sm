# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Create isolated run directories and write explicit, inspectable outputs."""

import csv
import hashlib
import json
import platform
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import matplotlib

from matplotlib.figure import Figure
import numpy as np
from PIL import __version__ as pillow_version

from .plotting import styled_plot
from . import __version__
from .model import InitialConfig, InitialState


# Output organisation and provenance.
def freeze_animation_layout(figure: Figure) -> None:
    """Settle automatic layout once, then keep plot and colourbar geometry fixed.

    Constrained layout can resize colourbars over successive draws even when
    their numeric limits are constant. Warm up before recording the first frame,
    then retain the resolved axes positions and disable automatic locators.
    """
    figure.canvas.draw()
    figure.canvas.draw()
    positions = [axes.get_position().frozen() for axes in figure.axes]
    figure.set_layout_engine("none")
    for axes, position in zip(figure.axes, positions):
        axes.set_axes_locator(None)
        axes.set_position(position)


def write_json(path: Path, value: object) -> None:
    """Write JSON atomically within the destination directory."""
    # Write to a neighbouring temporary file first, then replace the destination.
    # Readers see a complete JSON document rather than a partially written one.
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def create_run_directory(output_root: Path, mode: str = "initial") -> Path:
    """Create a unique child directory, preserving existing runs."""
    output_root = output_root.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    # UTC timestamps make runs easy to sort. A random suffix further reduces
    # name collisions; mkdir fails if a directory with that name already exists.
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    prefix = mode if mode in ("binding", "movement", "response", "full") else "initial-state"
    run_directory = output_root / f"{prefix}_{stamp}_{uuid4().hex[:8]}"
    run_directory.mkdir()
    return run_directory


def software_metadata() -> dict[str, object]:
    """Record software versions and hashes of the executing Python source files."""
    package_directory = Path(__file__).parent
    # A SHA-256 fingerprint changes when file bytes change, including comments.
    # This identifies the source used even before the package version changes.
    hashes = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(package_directory.glob("*.py"))
    }
    return {
        "program_version": __version__,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
        "pillow": pillow_version,
        "platform": platform.platform(),
        "random_bit_generator": "PCG64",
        "source_sha256": hashes,
    }


# Numerical outputs: retain the counts behind the plots for later inspection.
def save_state(directory: Path, state: InitialState) -> None:
    """Save numerical arrays without object pickles, plus a readable position table."""
    # NPZ stores named NumPy arrays together. Load with allow_pickle=False;
    # the state arrays contain only integers, so Python object pickles are unused.
    np.savez_compressed(
        directory / "initial_state.npz",
        occupancy=state.occupancy,
        free_signals=state.free_signals,
        free_polymers=state.free_polymers,
    )
    # One row per bacterium records the lower-left anchor of its 2 x 2 footprint.
    with (directory / "bacteria.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["identifier", "x", "y"])
        writer.writeheader()
        writer.writerows(asdict(bacterium) for bacterium in state.bacteria)


# Visual outputs: cell positions and free-entity counts share grid coordinates.
@styled_plot
def save_plot(directory: Path, state: InitialState, config: InitialConfig) -> None:
    """Draw initial counts and bacterial footprints; this is not a simulation result."""
    # Select a noninteractive backend before importing pyplot, allowing plots
    # to be saved on computers without a graphical desktop.
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import MaxNLocator

    size = config.grid_size
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), layout="constrained")
    try:
        axes[0].set_title(f"Bacterial positions (n = {len(state.bacteria)})")
        # Overlay the same 2 x 2 footprints on all panels. Outlines on heatmaps
        # show where bacteria sit without covering the underlying counts.
        for ax in axes:
            ax.set(xlim=(0, size), ylim=(0, size), xlabel="x (lattice units)")
            ax.set_ylabel("y (lattice units)")
            ax.set_aspect("equal")
            for bacterium in state.bacteria:
                ax.add_patch(
                    Rectangle(
                        (bacterium.x, bacterium.y), 2, 2,
                        facecolor="#0072B2" if ax is axes[0] else "none",
                        edgecolor="#333333" if ax is axes[0] else "#FFFFFF",
                        linewidth=0.6,
                    )
                )
        for ax, counts, title in (
            (axes[1], state.free_signals, "Free signal count"),
            (axes[2], state.free_polymers, "Free polymer count"),
        ):
            # origin='lower' aligns array row y=0 with the bottom of the plot.
            # Each heatmap has its own count scale; compare the colourbars.
            image = ax.imshow(
                counts, origin="lower", extent=(0, size, 0, size),
                interpolation="nearest", cmap="cividis", vmin=0,
                vmax=max(1, int(counts.max())), zorder=0,
            )
            ax.set_title(f"{title} (total = {int(counts.sum())})")
            colorbar = fig.colorbar(
                image, ax=ax, label="Count per lattice site", shrink=0.75
            )
            # Counts are discrete. An empty panel needs only a zero tick;
            # vmax remains at least 1 to keep a valid colour range.
            if counts.any():
                colorbar.locator = MaxNLocator(integer=True)
                colorbar.update_ticks()
            else:
                colorbar.set_ticks([0])
        fig.suptitle(f"{config.strain}: initial state only | seed {config.seed}")
        fig.savefig(directory / "initial_state.png", dpi=180)
        fig.savefig(directory / "initial_state.svg")
    finally:
        # Release the figure even if saving fails, especially across many runs.
        plt.close(fig)

