# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Delayed per-cell expression and predicted light from supplementary S21-S22.

RESP-01 through RESP-05 in model-notes.md specify the decay-sign interpretation,
piecewise solution, initial expression, light transformation, and aggregation.
This is a downstream response: it does not feed back into binding or movement.
"""

from dataclasses import dataclass
import math

import numpy as np
from numpy.typing import NDArray

from .binding import BindingState


MAX_SENSING = 60
LIGHT_SLOPE = 0.25


@dataclass(frozen=True)
class ResponseConfig:
    """Response controls; expression is a continuous proxy, not a molecule count."""

    k_syn: float = 0.05
    k_deg: float = 0.01
    response_initial: float = 0.0
    response_midpoint: float | None = None
    light_scale: float = 2500.0
    light_map_max: float | None = None

    def __post_init__(self) -> None:
        for name in ("k_syn", "k_deg", "response_initial", "response_midpoint",
                     "light_scale", "light_map_max"):
            value = getattr(self, name)
            if name in ("response_midpoint", "light_map_max") and value is None:
                continue
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be a finite non-negative number.")
        if self.response_initial > MAX_SENSING:
            raise ValueError("response_initial must be between 0 and 60.")
        if self.light_scale <= 0 or not math.isfinite(2.0 * self.light_scale):
            raise ValueError("light_scale and twice light_scale must be positive and finite.")
        if self.light_map_max is not None and self.light_map_max <= 0:
            raise ValueError("light_map_max must be positive when supplied.")

    def midpoint(self, strain: str) -> float:
        """Use the S22 strain midpoint unless an explicit user override is supplied."""
        if strain not in ("MM32", "BB170"):
            raise ValueError("strain must be MM32 or BB170.")
        return self.response_midpoint if self.response_midpoint is not None else (
            38.0 if strain == "MM32" else 15.0
        )


def advance_expression(
    expression: NDArray[np.float64], sensed: NDArray[np.int64],
    elapsed: NDArray[np.float64], config: ResponseConfig,
) -> NDArray[np.float64]:
    """Integrate the piecewise ODE exactly while receptor occupancy is constant.

    If sensing exceeds expression, E approaches S with rate k_syn. If E exceeds
    S, E decays exponentially until it reaches S, then remains there. Clamping
    at that crossing implements the branch change instead of overshooting it.
    """
    if np.any(~np.isfinite(elapsed)) or np.any(elapsed < 0):
        raise ValueError("Response elapsed times must be finite and non-negative.")
    if (np.any(~np.isfinite(expression)) or np.any(expression < 0)
            or np.any(expression > MAX_SENSING)
            or np.any(sensed < 0) or np.any(sensed > MAX_SENSING)):
        raise ValueError("Expression and receptor occupancy must lie between 0 and 60.")
    # expm1 preserves small increases when k_syn*elapsed is close to zero.
    rising = expression + (sensed - expression) * (-np.expm1(-config.k_syn * elapsed))
    falling = np.maximum(sensed, expression * np.exp(-config.k_deg * elapsed))
    return np.where(sensed >= expression, rising, falling)


def light_from_expression(
    expression: NDArray[np.float64], midpoint: float, scale: float,
) -> NDArray[np.float64]:
    """Evaluate the printed light curve stably, including its nonzero baseline.

    Light = 2*Scale / (1 + exp(0.25*(midpoint-E))). logaddexp avoids overflow
    without the cancellation that would occur in 1-tanh for very dim cells.
    """
    delta = LIGHT_SLOPE * (midpoint - expression)
    return (2.0 * scale) * np.exp(-np.logaddexp(0.0, delta))


class ResponseTracker:
    """Track expression lazily between actual bacterial signal-binding events.

    Each cell has its own last occupancy-change time. Observations evaluate the
    exact solution without modifying that clock. Extra samples/frames therefore
    do not change the result, and frequent unrelated signal hops do not require
    updating every cell's expression.
    """

    def __init__(self, state: BindingState, config: ResponseConfig, strain: str) -> None:
        self.config = config
        self.midpoint = config.midpoint(strain)
        self.identifiers = np.array([b.identifier for b in state.bacteria], dtype=np.int64)
        self.indices = {int(identifier): index for index, identifier in enumerate(self.identifiers)}
        self.sensed = self.sensing_from_state(state)
        self.expression = np.full(len(self.identifiers), config.response_initial, dtype=float)
        self.updated = np.zeros(len(self.identifiers))

    @staticmethod
    def sensing_from_state(state: BindingState) -> NDArray[np.int64]:
        """Sum receptor-bound signals over each cell's four current footprint sites."""
        return np.array([
            int(state.bacteria_signals[b.y:b.y + 2, b.x:b.x + 2].sum())
            for b in state.bacteria
        ], dtype=np.int64)

    def change_sensing(self, bacterium: int, change: int, time: float) -> None:
        """Advance under the old occupancy, then apply one binding/release change."""
        index = self.indices[bacterium]
        self.expression[index] = advance_expression(
            self.expression[index:index + 1], self.sensed[index:index + 1],
            np.array([time - self.updated[index]]), self.config,
        )[0]
        self.updated[index] = time
        self.sensed[index] += change
        if not 0 <= self.sensed[index] <= MAX_SENSING:
            raise ValueError("Response sensing exceeds receptor capacity.")

    def add_daughter(self, identifier: int, time: float, expression: float) -> None:
        """Start an unbound daughter at birth; the parent's response is retained.

        The source does not specify expression inheritance. POP-04 defaults to
        zero daughter expression, exposed separately from initial-cell expression.
        """
        if identifier in self.indices:
            raise ValueError("A daughter identifier is already present in the response tracker.")
        self.indices[identifier] = len(self.identifiers)
        self.identifiers = np.append(self.identifiers, identifier)
        self.sensed = np.append(self.sensed, 0)
        self.expression = np.append(self.expression, expression)
        self.updated = np.append(self.updated, time)

    def observe(self, time: float) -> dict[str, NDArray]:
        """Return the current response without changing internal integration anchors."""
        expression = advance_expression(self.expression, self.sensed,
                                        time - self.updated, self.config)
        light = light_from_expression(expression, self.midpoint, self.config.light_scale)
        return {"bacterium": self.identifiers.copy(), "sensed": self.sensed.copy(),
                "expression": expression, "luminescence": light}

    def validate(self, state: BindingState, time: float) -> None:
        """Check that receptor bookkeeping matches the spatial state and light is bounded."""
        if not np.array_equal(self.sensing_from_state(state), self.sensed):
            raise ValueError("Spatial receptor occupancy and response sensing disagree.")
        response = self.observe(time)
        light = response["luminescence"]
        if (np.any(~np.isfinite(light)) or np.any(light < 0)
                or np.any(light > 2 * self.config.light_scale)):
            raise ValueError("Predicted light is non-finite or outside its model bounds.")


def response_summary(response: dict[str, NDArray], time: float) -> dict[str, int | float]:
    """Report raw population sums and means, with zero means for an empty population."""
    count = len(response["bacterium"])
    return {
        "time": time, "bacteria": count,
        "sensed_total": int(response["sensed"].sum()),
        "sensed_mean": float(response["sensed"].mean()) if count else 0.0,
        "expression_total": float(response["expression"].sum()),
        "expression_mean": float(response["expression"].mean()) if count else 0.0,
        "luminescence_total": float(response["luminescence"].sum()),
        "luminescence_mean": float(response["luminescence"].mean()) if count else 0.0,
    }
