# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 21 October 2026

"""Reversible binding on a fixed lattice, using Gillespie's direct method.

Source: supplementary pages S20-S21 and S23. Geometry, local site allocation,
and kinetic choices are documented as BIND-01 through BIND-06 in model-notes.md.
This diagnostic stage has no movement, growth, signal production, or light output.
"""

from dataclasses import dataclass, field
import math
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from .model import Bacterium, InitialState


# Published totals are 60 signal sites and 120 polymer sites per bacterium.
# BIND-01 distributes them equally over its four occupied lattice sites.
BS_PER_SITE = 15
BP_PER_SITE = 30
PB_PER_POLYMER = 2
PS_PER_POLYMER = 5
CROWDING_THRESHOLD = 10


@dataclass(frozen=True)
class BindingConfig:
    """User controls for the fixed-position binding experiment (arbitrary units)."""

    duration: float = 100.0
    sample_interval: float = 5.0
    max_events: int = 500_000
    k_bs_on: float = 0.001
    k_bs_off: float = 0.002
    # The supplement does not identify these individual polymer rates. Require
    # explicit choices when polymers are present instead of inventing defaults.
    k_pb_on: float | None = None
    k_pb_off: float | None = None
    k_ps_on: float | None = None
    k_ps_off: float | None = None
    # S23's site-occlusion rule remains the default. None is an explicit
    # sensitivity intervention: retain all empty PS slots at any local density.
    ps_crowding_threshold: int | None = CROWDING_THRESHOLD

    def __post_init__(self) -> None:
        for name in ("duration", "sample_interval", "k_bs_on", "k_bs_off",
                     "k_pb_on", "k_pb_off", "k_ps_on", "k_ps_off"):
            value = getattr(self, name)
            if value is None and name.startswith(("k_pb_", "k_ps_")):
                continue
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be a finite non-negative number.")
        if self.duration <= 0 or self.sample_interval <= 0:
            raise ValueError("duration and sample_interval must be greater than zero.")
        if self.duration / self.sample_interval > 100_000:
            raise ValueError("Choose a sample_interval giving at most 100000 intervals.")
        if type(self.max_events) is not int or self.max_events < 1:
            raise ValueError("max_events must be a positive integer.")
        if self.ps_crowding_threshold is not None and (
            type(self.ps_crowding_threshold) is not int or self.ps_crowding_threshold < 1
        ):
            raise ValueError("ps_crowding_threshold must be a positive integer or null (off).")

    def validate_polymer_rates(self, polymers: int) -> None:
        """Require all four polymer rates, including explicit zeros to disable binding."""
        if polymers and any(getattr(self, name) is None for name in (
            "k_pb_on", "k_pb_off", "k_ps_on", "k_ps_off"
        )):
            raise ValueError(
                "Binding with polymers requires k_pb_on, k_pb_off, k_ps_on, and "
                "k_ps_off. Use examples/mm32_p3_binding.json for documented example "
                "choices; use explicit zeros for a disabled reaction."
            )

    def rate_metadata(self) -> dict[str, dict[str, float | str | None]]:
        """Record rates, ratios, and whether the individual kinetics are published."""
        records = {}
        for pair in ("bs", "pb", "ps"):
            on = getattr(self, f"k_{pair}_on")
            off = getattr(self, f"k_{pair}_off")
            # A zero reverse rate has no finite association ratio. JSON null
            # avoids writing Infinity, which is not valid standard JSON.
            ratio = on / off if on is not None and off is not None and off > 0 else None
            if ratio is not None and not math.isfinite(ratio):
                ratio = None
            origin = ("S21 published pair" if pair == "bs" and on == 0.001 and off == 0.002
                      else "user-selected kinetics" if on is not None or off is not None
                      else "unspecified; no polymers present")
            records[pair] = {"on": on, "off": off, "association_ratio": ratio,
                             "origin": origin}
        return records


@dataclass
class BindingState:
    """Counts plus explicit polymer identities and their two bacterial contacts.

    Fixed-position binding retains polymer anchors. Movement mode updates them
    using the carrier and tether rules in movement.py.
    A contact stores a flattened site index; -1 means an unoccupied polymer end.
    Bacterial identities at contacts are read from the occupancy grid.
    """

    occupancy: NDArray[np.int32]
    bacteria: tuple[Bacterium, ...]
    free_signals: NDArray[np.int64]
    bacteria_signals: NDArray[np.int64]
    bacteria_polymers: NDArray[np.int64]
    polymer_sites: NDArray[np.int64]
    polymer_signals: NDArray[np.int64]
    polymer_contacts: NDArray[np.int64]

    @classmethod
    def from_initial(cls, initial: InitialState) -> "BindingState":
        """Copy initial arrays; assign polymer IDs in row-major lattice order."""
        # Explicit polymer identity is needed to distinguish one polymer bridge
        # from two separate polymers, each attached to a different bacterium.
        sites = np.repeat(np.arange(initial.occupancy.size), initial.free_polymers.ravel())
        count = len(sites)
        return cls(
            occupancy=initial.occupancy.copy(), bacteria=initial.bacteria,
            free_signals=initial.free_signals.copy(),
            bacteria_signals=np.zeros_like(initial.free_signals),
            bacteria_polymers=np.zeros_like(initial.free_signals),
            polymer_sites=sites.astype(np.int64),
            polymer_signals=np.zeros(count, dtype=np.int64),
            polymer_contacts=np.full((count, PB_PER_POLYMER), -1, dtype=np.int64),
        )

    def validate(self, total_signals: int, total_polymers: int) -> None:
        """Check conservation, capacities, and agreement of both contact records."""
        shape = self.occupancy.shape
        # Recheck geometry after cluster movement as well as after chemistry.
        expected = np.full(shape, -1, dtype=np.int32)
        for bacterium in self.bacteria:
            x, y = bacterium.x, bacterium.y
            if not (0 <= x < shape[1] - 1 and 0 <= y < shape[0] - 1):
                raise ValueError("A bacterium lies outside the lattice.")
            footprint = expected[y:y + 2, x:x + 2]
            if np.any(footprint >= 0):
                raise ValueError("Bacteria overlap.")
            footprint[:] = bacterium.identifier
        if not np.array_equal(expected, self.occupancy):
            raise ValueError("Bacterial positions and occupancy disagree.")
        for counts, capacity in ((self.free_signals, None),
                                 (self.bacteria_signals, BS_PER_SITE),
                                 (self.bacteria_polymers, BP_PER_SITE)):
            if counts.shape != shape or np.any(counts < 0):
                raise ValueError("Invalid lattice counts or shape.")
            if capacity is not None:
                if np.any(counts > capacity) or np.any(counts[self.occupancy < 0]):
                    raise ValueError("Bacterial binding-site capacity was violated.")
        if (len(self.polymer_sites) != total_polymers
                or self.polymer_signals.shape != (total_polymers,)
                or self.polymer_contacts.shape != (total_polymers, PB_PER_POLYMER)):
            raise ValueError("Polymer identity or capacity records are inconsistent.")
        if np.any(self.polymer_sites < 0) or np.any(self.polymer_sites >= self.occupancy.size):
            raise ValueError("A polymer anchor lies outside the lattice.")
        if np.any(self.polymer_signals < 0) or np.any(self.polymer_signals > PS_PER_POLYMER):
            raise ValueError("Polymer signal-site capacity was violated.")
        signals = int(self.free_signals.sum() + self.bacteria_signals.sum()
                      + self.polymer_signals.sum())
        if signals != total_signals:
            raise ValueError("Signal conservation was violated.")
        contacts = self.polymer_contacts[self.polymer_contacts >= 0]
        if np.any(self.polymer_contacts < -1) or np.any(contacts >= self.occupancy.size):
            raise ValueError("A polymer contact lies outside the lattice.")
        if np.any(self.occupancy.ravel()[contacts] < 0):
            raise ValueError("A polymer is attached to an empty lattice site.")
        counted = np.bincount(contacts, minlength=self.occupancy.size).reshape(shape)
        if not np.array_equal(counted, self.bacteria_polymers):
            raise ValueError("Polymer and bacterial contact records disagree.")
        for first, second in self.polymer_contacts[np.all(self.polymer_contacts >= 0, axis=1)]:
            if self.occupancy.ravel()[first] == self.occupancy.ravel()[second]:
                raise ValueError("Both ends of a bridge must connect distinct bacteria.")


@dataclass(frozen=True)
class Reaction:
    """One channel: a reaction type, entity/site index, rate, and optional contact."""

    kind: str
    target: int
    propensity: float
    contact: int = -1


@dataclass
class BindingResult:
    """Final state, sampled counts, and the actual stochastic event history."""

    state: BindingState
    status: str
    time: float
    trajectory: list[dict[str, int | float]] = field(default_factory=list)
    events: list[dict[str, int | float | str]] = field(default_factory=list)


def neighbour_contacts(state: BindingState) -> dict[int, tuple[int, ...]]:
    """Find sites on other bacteria sharing an edge with each bacterial footprint.

    BIND-02: diagonal contact alone is insufficient. Each reachable target site
    appears once, even if it touches the host through more than one edge.
    """
    occupancy = state.occupancy
    size = occupancy.shape[0]
    neighbours: dict[int, tuple[int, ...]] = {}
    for bacterium in state.bacteria:
        sites: set[int] = set()
        for y in range(bacterium.y, bacterium.y + 2):
            for x in range(bacterium.x, bacterium.x + 2):
                for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                    ny, nx = y + dy, x + dx
                    if (0 <= ny < size and 0 <= nx < size
                            and occupancy[ny, nx] >= 0
                            and occupancy[ny, nx] != bacterium.identifier):
                        sites.add(ny * size + nx)
        neighbours[bacterium.identifier] = tuple(sorted(sites))
    return neighbours


def bacterial_reactions(state: BindingState, config: BindingConfig, site: int) -> list[Reaction]:
    """Return the two local BS channels; used by reference and cached solvers."""
    channels: list[Reaction] = []
    if state.occupancy.ravel()[site] < 0:
        return channels
    bound = int(state.bacteria_signals.ravel()[site])
    for kind, rate in (
        ("bs_bind", config.k_bs_on * (BS_PER_SITE - bound)
         * int(state.free_signals.ravel()[site])),
        ("bs_unbind", config.k_bs_off * bound),
    ):
        if not math.isfinite(rate) or rate < 0:
            raise ValueError("A propensity is negative or non-finite.")
        if rate > 0:
            channels.append(Reaction(kind, site, rate))
    return channels


def polymer_reactions(
    state: BindingState, config: BindingConfig, neighbours: dict[int, tuple[int, ...]],
    polymer: int, local_count: int,
) -> list[Reaction]:
    """Return PS/PB channels for one polymer, retaining its contact identity."""
    channels: list[Reaction] = []
    site = int(state.polymer_sites[polymer])
    occupied = state.occupancy.ravel()

    def add(kind: str, target: int, rate: float, contact: int = -1) -> None:
        if not math.isfinite(rate) or rate < 0:
            raise ValueError("A propensity is negative or non-finite; check rate magnitudes.")
        if rate > 0:
            channels.append(Reaction(kind, target, rate, contact))

    bound_signals = int(state.polymer_signals[polymer])
    threshold = config.ps_crowding_threshold
    accessible = (PS_PER_POLYMER - bound_signals
                  if threshold is None or local_count < threshold else 0)
    add("ps_bind", polymer, (config.k_ps_on or 0.0) * accessible
        * int(state.free_signals.ravel()[site]))
    add("ps_unbind", polymer, (config.k_ps_off or 0.0) * bound_signals)
    contacts = state.polymer_contacts[polymer]
    used = np.flatnonzero(contacts >= 0)
    for end in used:
        add("pb_unbind", polymer, config.k_pb_off or 0.0, int(end))
    if len(used) == 0:
        targets = (site,) if occupied[site] >= 0 else ()
    elif len(used) == 1:
        host = int(occupied[contacts[used[0]]])
        targets = neighbours[host]
    else:
        targets = ()
    for target in targets:
        add("pb_bind", polymer, (config.k_pb_on or 0.0) * (PB_PER_POLYMER - len(used))
            * (BP_PER_SITE - int(state.bacteria_polymers.ravel()[target])), target)
    return channels


def reaction_channels(
    state: BindingState, config: BindingConfig,
    neighbours: dict[int, tuple[int, ...]],
) -> list[Reaction]:
    """Construct all chemical channels for the readable reference solver."""
    channels = []
    for site in np.flatnonzero(state.occupancy.ravel() >= 0):
        channels.extend(bacterial_reactions(state, config, int(site)))
    local = np.bincount(state.polymer_sites, minlength=state.occupancy.size)
    for polymer, site in enumerate(state.polymer_sites):
        channels.extend(polymer_reactions(state, config, neighbours, polymer, int(local[site])))
    return channels


def apply_reaction(state: BindingState, reaction: Reaction) -> dict[str, int | str]:
    """Apply one valid channel and return identifiers for the event log.

    Binding moves a signal out of the free pool; unbinding returns it to the
    corresponding site. Neither operation creates or destroys a signal.
    """
    target = reaction.target
    kind = reaction.kind
    polymer = -1
    occupied = state.occupancy.ravel()
    if kind in ("bs_bind", "bs_unbind"):
        site = target
        change = 1 if kind == "bs_bind" else -1
        state.free_signals.ravel()[site] -= change
        state.bacteria_signals.ravel()[site] += change
    elif kind in ("ps_bind", "ps_unbind"):
        polymer = target
        site = int(state.polymer_sites[polymer])
        change = 1 if kind == "ps_bind" else -1
        state.free_signals.ravel()[site] -= change
        state.polymer_signals[polymer] += change
    elif kind == "pb_bind":
        polymer = target
        site = reaction.contact
        end = int(np.flatnonzero(state.polymer_contacts[polymer] < 0)[0])
        state.polymer_contacts[polymer, end] = site
        state.bacteria_polymers.ravel()[site] += 1
    elif kind == "pb_unbind":
        polymer = target
        site = int(state.polymer_contacts[polymer, reaction.contact])
        state.polymer_contacts[polymer, reaction.contact] = -1
        state.bacteria_polymers.ravel()[site] -= 1
    else:
        raise ValueError(f"Unknown reaction: {kind}")
    y, x = divmod(site, state.occupancy.shape[0])
    return {"reaction": kind, "x": x, "y": y, "polymer": polymer,
            "bacterium": int(occupied[site]) if kind.startswith(("bs_", "pb_")) else -1}


def sample_counts(state: BindingState, time: float) -> dict[str, int | float]:
    """Summarise signal pools and the polymer-mediated bacterial contact graph."""
    attachments = np.count_nonzero(state.polymer_contacts >= 0, axis=1)
    # Union-find identifies connected bacteria. Parallel bridges between the
    # same pair still count as two bridges, but do not create extra bacteria.
    parents = {b.identifier: b.identifier for b in state.bacteria}

    def find(identifier: int) -> int:
        while parents[identifier] != identifier:
            identifier = parents[identifier]
        return identifier

    for first, second in state.polymer_contacts[attachments == 2]:
        a, b = int(state.occupancy.ravel()[first]), int(state.occupancy.ravel()[second])
        parents[find(a)] = find(b)
    sizes: dict[int, int] = {}
    for identifier in parents:
        root = find(identifier)
        sizes[root] = sizes.get(root, 0) + 1
    free = int(state.free_signals.sum())
    bs = int(state.bacteria_signals.sum())
    ps = int(state.polymer_signals.sum())
    return {
        "time": time, "free_signals": free, "bacteria_bound_signals": bs,
        "polymer_bound_signals": ps, "total_signals": free + bs + ps,
        "free_polymers": int(np.count_nonzero(attachments == 0)),
        "attached_polymers": int(np.count_nonzero(attachments > 0)),
        "pb_bonds": int(attachments.sum()), "bridges": int(np.count_nonzero(attachments == 2)),
        "largest_cluster": max(sizes.values(), default=0),
    }


def simulate_binding(initial: InitialState, config: BindingConfig, seed: int,
                     progress: Callable[[float, int], None] | None = None) -> BindingResult:
    """Run exact direct SSA for this documented fixed-position reaction system.

    Observation times do not alter event timing. If the event cap is reached,
    return an explicitly incomplete trajectory ending at the actual stop time.
    An optional read-only callback observes actual time and event count at start,
    every 1000 events, and stop; it introduces no observation times or RNG draws.
    """
    if type(seed) is not int or not 0 <= seed < 2**64:
        raise ValueError("seed must be an integer between 0 and 2**64 - 1.")
    total_signals = int(initial.free_signals.sum())
    total_polymers = int(initial.free_polymers.sum())
    config.validate_polymer_rates(total_polymers)
    state = BindingState.from_initial(initial)
    state.validate(total_signals, total_polymers)
    # BIND-06: a distinct deterministic stream avoids restarting the placement
    # generator's sequence. No global random state is read or modified.
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([seed, 1])))
    neighbours = neighbour_contacts(state)
    result = BindingResult(state, "complete", 0.0, [sample_counts(state, 0.0)])
    next_index = 1
    if progress is not None:
        progress(0.0, 0)

    def observations_before(limit: float, inclusive: bool) -> None:
        nonlocal next_index
        while True:
            time = next_index * config.sample_interval
            if time > config.duration or (time > limit if inclusive else time >= limit):
                break
            result.trajectory.append(sample_counts(state, time))
            next_index += 1

    while result.time < config.duration:
        channels = reaction_channels(state, config, neighbours)
        total = math.fsum(channel.propensity for channel in channels)
        if not math.isfinite(total):
            raise ValueError("Total propensity is non-finite; reduce rate magnitudes.")
        next_time = (result.time + float(rng.exponential(1.0 / total))
                     if total > 0 else math.inf)
        if next_time > config.duration:
            observations_before(config.duration, inclusive=True)
            result.time = config.duration
            break
        if next_time <= result.time:
            raise ValueError("Event time cannot advance at this numerical precision.")
        observations_before(next_time, inclusive=False)
        cumulative = np.cumsum([channel.propensity for channel in channels])
        # Use the same cumulative sum for the draw and search, avoiding a
        # floating-point discrepancy between different ways of adding rates.
        choice = float(rng.random()) * float(cumulative[-1])
        selected = min(int(np.searchsorted(cumulative, choice, side="right")), len(channels) - 1)
        record = apply_reaction(state, channels[selected])
        result.time = next_time
        result.events.append({"event": len(result.events) + 1, "time": next_time, **record})
        if progress is not None and len(result.events) % 1000 == 0:
            progress(result.time, len(result.events))
        state.validate(total_signals, total_polymers)
        observations_before(next_time, inclusive=True)
        if len(result.events) >= config.max_events and result.time < config.duration:
            result.status = "event_limit"
            break
    if result.trajectory[-1]["time"] != result.time:
        result.trajectory.append(sample_counts(state, result.time))
    state.validate(total_signals, total_polymers)
    if progress is not None:
        progress(result.time, len(result.events))
    return result
