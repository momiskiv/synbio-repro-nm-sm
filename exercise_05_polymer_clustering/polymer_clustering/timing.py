# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Measure real execution time independently of the model's arbitrary time units."""

from contextlib import contextmanager
from time import perf_counter, process_time
from typing import Iterator


class RunTimer:
    """Accumulate non-overlapping run phases using a monotonic elapsed-time clock.

    Wall time measures how long the user waits, including time the process is
    not scheduled. CPU time counts work done by this process, rather than model
    time. Neither measurement makes random draws or enters the event equations.
    """

    def __init__(self) -> None:
        self.wall_started = perf_counter()
        self.cpu_started = process_time()
        self.phases: dict[str, float] = {}
        self.active_phase: str | None = None

    @contextmanager
    def phase(self, name: str) -> Iterator[None]:
        """Include partially completed work even when a phase raises an exception."""
        if self.active_phase is not None:
            raise ValueError("Timing phases must not overlap.")
        self.active_phase = name
        started = perf_counter()
        try:
            yield
        finally:
            self.phases[name] = self.phases.get(name, 0.0) + perf_counter() - started
            self.active_phase = None

    def snapshot(self, events: int | None = None) -> dict[str, object]:
        """Record seconds up to this snapshot; final metadata writing follows it."""
        total = perf_counter() - self.wall_started
        simulation = self.phases.get("simulation", 0.0)
        return {
            "clock": "time.perf_counter",
            "total_wall_seconds": total,
            "cpu_seconds": process_time() - self.cpu_started,
            "initialisation_wall_seconds": self.phases.get("initialisation", 0.0),
            "simulation_wall_seconds": simulation,
            "output_wall_seconds": self.phases.get("output", 0.0),
            "other_wall_seconds": max(0.0, total - sum(self.phases.values())),
            "simulation_events_per_second": (
                events / simulation if events is not None and simulation > 0 else None
            ),
            "scope": (
                "after argument parsing to final metadata snapshot; "
                "excludes imports, final metadata write, and terminal printing"
            ),
            "cpu_scope": "current process only; excludes any external worker processes",
        }


def format_elapsed(seconds: float) -> str:
    """Display real seconds and, for longer runs, an hours:minutes:seconds reading."""
    text = f"{seconds:.3f} s"
    if seconds >= 60:
        milliseconds = round(seconds * 1000)
        hours, remainder = divmod(milliseconds, 3600000)
        minutes, remainder = divmod(remainder, 60000)
        text += f" ({hours:02d}:{minutes:02d}:{remainder / 1000:06.3f})"
    return text


def timing_summary(timing: dict[str, object]) -> str:
    """Keep terminal reporting short while preserving detailed seconds in run.json."""
    return (
        f"Runtime: total {format_elapsed(float(timing['total_wall_seconds']))}; "
        f"simulation {float(timing['simulation_wall_seconds']):.3f} s; "
        f"outputs {float(timing['output_wall_seconds']):.3f} s."
    )
