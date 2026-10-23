# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Initial spatial state; no reactions, movement, growth, or light calculation yet.

Paper-derived rules and reconstruction choices are listed in docs/model-notes.md.
Array indexing is [y, x]; coordinates denote a bacterium's lower-left lattice site.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


# Configuration: values describing the requested world before it is created.
# A dataclass supplies an initialiser; frozen=True prevents reassignment of fields.
@dataclass(frozen=True)
class InitialConfig:
    """User-controlled initial conditions, validated before a run is created."""

    grid_size: int = 50
    bacteria: int = 25
    polymers: int = 0
    strain: str = "MM32"
    seed: int = 2026
    rng_scheme: str = "independent_v1"

    def __post_init__(self) -> None:
        # Validate immediately after construction, before generating any outputs.
        # Python treats bool as a subclass of int; an exact type check rejects it.
        for name in ("grid_size", "bacteria", "polymers", "seed"):
            if type(getattr(self, name)) is not int:
                raise ValueError(f"{name} must be an integer (not a boolean).")
        if self.grid_size < 2:
            raise ValueError("grid_size must be at least 2 for a 2 x 2 bacterium.")
        if not 0 <= self.bacteria <= 80:
            raise ValueError("bacteria must be between 0 and the published model cap of 80.")
        # Each bacterium occupies 2 x 2 sites. This is the geometric capacity;
        # sequential random placement can jam before that capacity is reached.
        capacity = (self.grid_size // 2) ** 2
        if self.bacteria > capacity:
            raise ValueError(f"This grid can hold at most {capacity} non-overlapping bacteria.")
        if self.polymers < 0:
            raise ValueError("polymers must be non-negative.")
        if not 0 <= self.seed < 2**64:
            raise ValueError("seed must be between 0 and 2**64 - 1.")
        if self.strain not in ("MM32", "BB170"):
            raise ValueError("strain must be MM32 or BB170.")
        if self.rng_scheme not in ("independent_v1", "legacy_v02"):
            raise ValueError("rng_scheme must be independent_v1 or legacy_v02.")

    @property
    def signal_count(self) -> int:
        """Initial total signal count stated on supplementary page S21."""
        # MM32 receives signal initially; BB170 starts without it. Signal
        # production by BB170 belongs to the future dynamics implementation.
        return 10_000 if self.strain == "MM32" else 0


# Spatial state: bacterial coordinates plus site-by-site counts of small entities.
@dataclass(frozen=True)
class Bacterium:
    """An identifier and the lower-left site of a fixed 2 x 2 footprint."""

    identifier: int
    x: int
    y: int


@dataclass
class InitialState:
    """Grid arrays and bacteria; -1 denotes an unoccupied bacterial site."""

    # Arrays use [row, column] = [y, x]. Occupancy stores bacterial identifiers;
    # the two count arrays allow multiple small entities to share any site.
    occupancy: NDArray[np.int32]
    bacteria: tuple[Bacterium, ...]
    free_signals: NDArray[np.int64]
    free_polymers: NDArray[np.int64]


def initialise(config: InitialConfig) -> InitialState:
    """Place bacteria sequentially at random valid anchors, then small entities.

    Reconstruction choices: closed edges, integer anchors, initially free signals
    and polymers, and independent uniform allocation of small entities to sites.
    This is not a uniform sample over all possible bacterial configurations.
    Dense configurations may jam; failure is explicit rather than changing rules.
    """
    # A local generator avoids changing NumPy's global random state. Repeating
    # the settings and seed in the same software environment repeats the draws.
    rng = np.random.Generator(np.random.PCG64(config.seed))
    # INIT-05: separate streams keep bacteria and signals identical when only
    # polymer count changes. Legacy mode retains the v0.2 draw order for replay.
    if config.rng_scheme == "independent_v1":
        placement_rng, polymer_rng, signal_rng = (
            np.random.Generator(np.random.PCG64(np.random.SeedSequence([config.seed, 0, i])))
            for i in range(3)
        )
    else:
        placement_rng = polymer_rng = signal_rng = rng
    size = config.grid_size
    occupancy = np.full((size, size), -1, dtype=np.int32)
    bacteria: list[Bacterium] = []
    # Valid lower-left anchors run from 0 to size - 2, keeping every 2 x 2
    # footprint inside the grid. The edges do not wrap around.
    anchors_per_axis = size - 1

    if config.bacteria:
        # Visit each possible anchor once, in a random order. divmod converts
        # the flattened anchor index into its row (y) and column (x).
        for anchor in placement_rng.permutation(anchors_per_axis**2):
            y, x = divmod(int(anchor), anchors_per_axis)
            # Slicing stops before y + 2 and x + 2, giving exactly four sites.
            # This slice is a view: writing to it also changes occupancy.
            footprint = occupancy[y : y + 2, x : x + 2]
            if np.any(footprint != -1):
                continue
            # Reject overlapping placements; accepted bacteria get IDs 0, 1, ...
            identifier = len(bacteria)
            footprint[:] = identifier
            bacteria.append(Bacterium(identifier, x, y))
            if len(bacteria) == config.bacteria:
                break
    # A crowded random arrangement can block all remaining anchors. Report the
    # failure explicitly rather than returning fewer bacteria than requested.
    if len(bacteria) != config.bacteria:
        raise ValueError(
            "Sequential random placement jammed before placing every bacterium. "
            "Use a different seed, fewer bacteria, or a larger grid."
        )

    # Our documented placement choice: each small entity independently has an
    # equal chance of occupying any site, including sites containing bacteria.
    # A multinomial draw allocates the exact requested total among those sites;
    # reshape restores the 2D grid without changing any counts.
    probabilities = np.full(size * size, 1.0 / (size * size))
    free_polymers = polymer_rng.multinomial(config.polymers, probabilities).reshape(size, size)
    free_signals = signal_rng.multinomial(config.signal_count, probabilities).reshape(size, size)
    # Store positions as a tuple and counts as explicit integer arrays. These
    # arrays describe the starting state; no binding or movement has occurred.
    return InitialState(
        occupancy=occupancy,
        bacteria=tuple(bacteria),
        free_signals=free_signals.astype(np.int64),
        free_polymers=free_polymers.astype(np.int64),
    )

