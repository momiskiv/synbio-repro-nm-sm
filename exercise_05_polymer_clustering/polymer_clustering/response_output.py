# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Save receptor, expression, and predicted-light outputs without conflating them."""

import csv
from pathlib import Path

import matplotlib

import numpy as np
from PIL import Image

from .plotting import styled_plot
from .movement import MovementResult
from .output import freeze_animation_layout


def save_response(
    directory: Path, result: MovementResult, light_map_max: float | None = None,
) -> None:
    """Save population/per-cell tables, final response, plots, and a paired map GIF."""
    if result.response_final is None:
        raise ValueError("Response outputs require a response-enabled simulation.")
    np.savez_compressed(directory / "response_state.npz", time=result.time,
                        **result.response_final)
    np.savez_compressed(
        directory / "response_frames.npz",
        times=np.array([frame.time for frame in result.frames]),
        bacterium=result.response_final["bacterium"],
        bacteria_counts=np.array([len(frame.bacteria) for frame in result.frames]),
        present=np.array([[index < len(frame.bacteria) for index in range(len(result.state.bacteria))]
                          for frame in result.frames], dtype=bool),
        sensed=padded_response(result, "sensed"),
        expression=padded_response(result, "expression"),
        luminescence=padded_response(result, "luminescence"),
    )
    for filename, rows, fields in (
        ("response_trajectory.csv", result.response_trajectory,
         list(result.response_trajectory[0])),
        ("cell_response.csv", result.cell_response,
         ["time", "bacterium", "x", "y", "sensed", "expression", "luminescence"]),
    ):
        with (directory / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    save_response_plot(directory, result)
    save_response_maps(directory, result, light_map_max)


def padded_response(result: MovementResult, name: str) -> np.ndarray:
    """Store absent cells as NaN, distinguishing them from real zero responses."""
    values = np.full((len(result.frames), len(result.state.bacteria)), np.nan)
    for index, frame in enumerate(result.frames):
        values[index, :len(frame.bacteria)] = getattr(frame, name)
    return values


@styled_plot
def save_response_plot(directory: Path, result: MovementResult) -> None:
    """Compare sampled mean receptor occupancy, expression proxy, and predicted light."""
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    rows = result.response_trajectory
    times = [row["time"] for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), layout="constrained")
    suffix = " | INCOMPLETE" if result.status != "complete" else ""
    if result.scenario_note:
        suffix += "\n" + result.scenario_note
    try:
        axes[0].step(times, [row["sensed_mean"] for row in rows], where="post",
                     label="Receptor occupancy", color="#0072B2")
        axes[0].plot(times, [row["expression_mean"] for row in rows],
                     label="Expression proxy", color="#D55E00")
        axes[0].set(title="Sensing and delayed expression", ylabel="Mean per-cell model state")
        axes[0].legend()
        axes[1].plot(times, [row["luminescence_mean"] for row in rows], color="#009E73")
        axes[1].set(title="Predicted light", ylabel="Mean light per cell (model a.u.)")
        for ax in axes:
            ax.set_xlabel("Time (arbitrary units)")
            ax.set_ylim(bottom=0)
            ax.grid(alpha=0.2)
        label = "Population-coupled response" if result.population_trajectory else (
            "Downstream response | no growth or signal production")
        fig.suptitle(label + suffix)
        fig.savefig(directory / "response.png", dpi=180)
        fig.savefig(directory / "response.svg")
    finally:
        plt.close(fig)


@styled_plot
def save_response_maps(
    directory: Path, result: MovementResult, light_map_max: float | None = None,
) -> None:
    """Pair free-signal maps with cell-specific light; keep scales fixed across frames."""
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.patches import Rectangle

    size = result.state.occupancy.shape[0]
    max_signal = max(1, max(int(frame.free_signals.max()) for frame in result.frames))
    max_light = light_map_max if light_map_max is not None else max(
        1.0, max(float(frame.luminescence.max(initial=0)) for frame in result.frames)
    )
    norm = Normalize(vmin=0, vmax=max_light)
    colours = plt.get_cmap("magma")
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), layout="constrained")
    pictures: list[Image.Image] = []
    suffix = " | INCOMPLETE" if result.status != "complete" else ""
    if result.scenario_note:
        suffix += "\n" + result.scenario_note
    try:
        for ax in axes:
            ax.set(xlim=(0, size), ylim=(0, size), xlabel="x (lattice units)",
                   ylabel="y (lattice units)")
            ax.set_aspect("equal")
        image = axes[0].imshow(result.frames[0].free_signals, origin="lower",
                               extent=(0, size, 0, size), interpolation="nearest",
                               cmap="cividis", vmin=0, vmax=max_signal)
        axes[0].set_title("Free signals and bacterial outlines")
        axes[1].set_title("Predicted light per bacterium")
        fig.colorbar(image, ax=axes[0], label="Free signals per site", shrink=0.75)
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=colours), ax=axes[1],
                     label="Predicted light per cell (model a.u.)", shrink=0.75)
        outlines, cells = [], []
        for bacterium in result.state.bacteria:
            outline = Rectangle((bacterium.x, bacterium.y), 2, 2, facecolor="none",
                                edgecolor="white", linewidth=0.8)
            cell = Rectangle((bacterium.x, bacterium.y), 2, 2,
                             edgecolor="#444444", linewidth=0.5)
            axes[0].add_patch(outline)
            axes[1].add_patch(cell)
            outline.set_visible(False)
            cell.set_visible(False)
            outlines.append(outline)
            cells.append(cell)
        fig.suptitle("Model response | t = 0 a.u." + suffix)
        freeze_animation_layout(fig)
        for frame in result.frames:
            image.set_data(frame.free_signals)
            for index in range(len(cells)):
                outlines[index].set_visible(index < len(frame.bacteria))
                cells[index].set_visible(index < len(frame.bacteria))
            for index, bacterium in enumerate(frame.bacteria):
                outlines[index].set_xy((bacterium.x, bacterium.y))
                cells[index].set_xy((bacterium.x, bacterium.y))
                cells[index].set_facecolor(colours(norm(frame.luminescence[index])))
            fig.suptitle(f"Model response | t = {frame.time:.3g} a.u." + suffix)
            fig.canvas.draw()
            pixels = np.asarray(fig.canvas.buffer_rgba()).copy()
            pictures.append(Image.fromarray(pixels).convert("RGB"))
        # Static views show the final observed state; GIF titles give time. The
        # light colourbar uses this run's peak cell value unless explicitly fixed.
        fig.savefig(directory / "response_spatial.png", dpi=180)
        fig.savefig(directory / "response_spatial.svg")
        pictures[0].save(directory / "response.gif", save_all=True,
                         append_images=pictures[1:], duration=200, loop=0)
    finally:
        plt.close(fig)
        for picture in pictures:
            picture.close()
