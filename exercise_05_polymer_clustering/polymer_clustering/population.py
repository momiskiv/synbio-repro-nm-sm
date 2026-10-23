# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Age-dependent division and BB170 signal production (supplement S19, S23).

POP-01 through POP-06 in model-notes.md distinguish source rules from choices.
Division times are sampled lifetimes, not exponential SSA reaction channels.
"""

from dataclasses import dataclass
import heapq
import math

import numpy as np

from .binding import BindingState
from .model import Bacterium


@dataclass(frozen=True)
class PopulationConfig:
    """Population controls, expressed in the model's arbitrary time units."""

    growth: bool = True
    doubling_time: float = 15000.0
    division_width: float = 20000.0
    population_cap: int = 80
    division_retry: float = 100.0
    production_rate: float | None = None
    daughter_expression: float = 0.0

    def __post_init__(self) -> None:
        if type(self.growth) is not bool:
            raise ValueError("growth must be true or false.")
        for name in ("doubling_time", "division_width", "division_retry",
                     "production_rate", "daughter_expression"):
            value = getattr(self, name)
            if name == "production_rate" and value is None:
                continue
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be a finite non-negative number.")
        if self.doubling_time <= 0 or self.doubling_time < self.division_width / 2:
            raise ValueError("doubling_time must be positive and at least half division_width.")
        if not math.isfinite(self.doubling_time + self.division_width / 2):
            raise ValueError("The division distribution upper bound must be finite.")
        if self.division_retry <= 0:
            raise ValueError("division_retry must be positive.")
        if type(self.population_cap) is not int or not 0 <= self.population_cap <= 80:
            raise ValueError("population_cap must be an integer between 0 and 80.")
        if self.daughter_expression > 60:
            raise ValueError("daughter_expression must be between 0 and 60.")

    def signal_rate(self, strain: str) -> float:
        """Default to 0.005 per BB170 cell, zero for MM32; overrides are experiments."""
        if strain not in ("MM32", "BB170"):
            raise ValueError("strain must be MM32 or BB170.")
        return self.production_rate if self.production_rate is not None else (
            0.005 if strain == "BB170" else 0.0
        )

    def lifetime(self, rng: np.random.Generator) -> float:
        """Interpret width as full support, symmetric about the stated mean.

        Width zero is a useful deterministic diagnostic. It avoids passing equal
        bounds to numpy's triangular sampler, which requires distinct bounds.
        """
        half = self.division_width / 2
        if half == 0:
            return self.doubling_time
        # S42 uses mean 10000; retaining S23's width 20000 puts the lower
        # support endpoint at zero. A continuous triangle has no atom there.
        # Redraw a finite-precision zero so appointments still advance in time.
        # Positive-support configurations keep exactly one draw, as before.
        value = 0.0
        while value <= 0.0:
            value = float(rng.triangular(self.doubling_time - half, self.doubling_time,
                                         self.doubling_time + half))
        return value


def nearest_empty_anchor(
    state: BindingState, parent: Bacterium, rng: np.random.Generator,
) -> tuple[int, int] | None:
    """Choose uniformly among nearest valid daughter anchors, using squared distance.

    Only bacterial occupancy blocks placement. Signals and polymers may already
    occupy an empty footprint and remain there when a daughter is introduced.
    """
    size = state.occupancy.shape[0]
    candidates: list[tuple[int, int]] = []
    nearest = math.inf
    for y in range(size - 1):
        for x in range(size - 1):
            if np.any(state.occupancy[y:y + 2, x:x + 2] >= 0):
                continue
            distance = (x - parent.x)**2 + (y - parent.y)**2
            if distance < nearest:
                nearest, candidates = distance, [(x, y)]
            elif distance == nearest:
                candidates.append((x, y))
    return candidates[int(rng.integers(len(candidates)))] if candidates else None


class DivisionSchedule:
    """Keep one due time per cell; reset both lifetimes only after successful division.

    All initial cells start at age zero (an explicit assumption). A heap permits
    deterministic division appointments to interrupt the stochastic SSA clock.
    Its identifier tie-break is reproducible when diagnostic lifetimes coincide.
    """

    def __init__(self, state: BindingState, config: PopulationConfig,
                 rng: np.random.Generator) -> None:
        if len(state.bacteria) > config.population_cap:
            raise ValueError("Initial bacteria exceed population_cap.")
        self.config, self.rng = config, rng
        self.queue: list[tuple[float, int]] = []
        self.lineage = [{"bacterium": b.identifier, "parent": -1, "birth_time": 0.0,
                         "generation": 0} for b in state.bacteria]
        self.generations = {b.identifier: 0 for b in state.bacteria}
        if config.growth and len(state.bacteria) < config.population_cap:
            for bacterium in state.bacteria:
                self.schedule(bacterium.identifier, 0.0)

    def schedule(self, identifier: int, time: float) -> None:
        """Draw a fresh lifetime, with an explicit failure if float time cannot advance."""
        self.push(identifier, time, self.config.lifetime(self.rng))

    def push(self, identifier: int, time: float, interval: float) -> None:
        due = time + interval
        if not math.isfinite(due) or due <= time:
            raise ArithmeticError("Division time cannot advance at this numerical precision.")
        heapq.heappush(self.queue, (due, identifier))

    @property
    def next_time(self) -> float:
        return self.queue[0][0] if self.queue else math.inf

    def divide(self, state: BindingState, time: float,
               placement_rng: np.random.Generator) -> dict[str, int | str]:
        """Place an unbound daughter without redistributing the parent's attachments."""
        due, identifier = heapq.heappop(self.queue)
        if due != time:
            raise ValueError("Division appointment and simulation clock disagree.")
        parent = next(b for b in state.bacteria if b.identifier == identifier)
        anchor = nearest_empty_anchor(state, parent, placement_rng)
        record: dict[str, int | str] = {
            "reaction": "division_blocked", "bacterium": identifier, "polymer": -1,
            "x": parent.x, "y": parent.y,
        }
        if anchor is None:
            # The supplement gives no failed-placement rule. Keep the overdue
            # parent ready and retry after a recorded positive interval (POP-03).
            self.push(identifier, time, self.config.division_retry)
            return record
        daughter = Bacterium(max(b.identifier for b in state.bacteria) + 1, *anchor)
        state.bacteria += (daughter,)
        state.occupancy[daughter.y:daughter.y + 2, daughter.x:daughter.x + 2] = daughter.identifier
        generation = self.generations[identifier] + 1
        self.generations[daughter.identifier] = generation
        self.lineage.append({"bacterium": daughter.identifier, "parent": identifier,
                             "birth_time": time, "generation": generation})
        if len(state.bacteria) >= self.config.population_cap:
            self.queue.clear()
        else:
            self.schedule(identifier, time)
            self.schedule(daughter.identifier, time)
        record.update(reaction="division", daughter=daughter.identifier,
                      to_x=daughter.x, to_y=daughter.y)
        return record


def produce_signal(state: BindingState, rng: np.random.Generator) -> dict[str, int | str]:
    """Choose an equal-rate producer, then one of its four sites uniformly."""
    bacterium = state.bacteria[int(rng.integers(len(state.bacteria)))]
    offset = int(rng.integers(4))
    x, y = bacterium.x + offset % 2, bacterium.y + offset // 2
    state.free_signals[y, x] += 1
    return {"reaction": "signal_production", "bacterium": bacterium.identifier,
            "polymer": -1, "x": x, "y": y}
