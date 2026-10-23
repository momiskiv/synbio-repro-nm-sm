# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Coupled binding and spatial transport, with explicit reconstruction choices.

Source rules: supplementary S19-S23. MOVE-01 through MOVE-07 in model-notes.md
specify obstruction, carrier transport, cluster displacement, and solver choices.
Optional population mode adds scheduled division and Poisson signal production. A downstream
response tracker evaluates expression and light without changing event rates.
"""

from collections import Counter
from dataclasses import dataclass, field
import math
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from .binding import (
    BS_PER_SITE, BindingConfig, BindingResult, BindingState, Reaction,
    apply_reaction, bacterial_reactions, neighbour_contacts, polymer_reactions, sample_counts,
)
from .model import Bacterium, InitialState
from .population import DivisionSchedule, PopulationConfig, produce_signal
from .rates import RateTree
from .response import ResponseConfig, ResponseTracker, response_summary


# The order is shared by neighbour tables, event records, and cluster updates.
DIRECTIONS = ((-1, 0), (1, 0), (0, -1), (0, 1))  # (dy, dx)


@dataclass(frozen=True)
class MovementConfig:
    """Movement rates and diagnostic controls; all rates use arbitrary units."""

    d_s: float = 0.05
    d_p: float = 0.01
    d_b: float = 0.0005
    obstruction_factor: float = 0.1
    crowding_lambda: float = 8.0
    advection: bool = True
    validate_every: int = 1000
    spatial_frames: int = 21
    event_log: str = "compact"

    def __post_init__(self) -> None:
        for name in ("d_s", "d_p", "d_b", "obstruction_factor", "crowding_lambda"):
            value = getattr(self, name)
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be a finite non-negative number.")
        if self.obstruction_factor > 1:
            raise ValueError("obstruction_factor must lie between 0 and 1.")
        if type(self.advection) is not bool:
            raise ValueError("advection must be true or false.")
        if type(self.validate_every) is not int or self.validate_every < 1:
            raise ValueError("validate_every must be a positive integer.")
        if type(self.spatial_frames) is not int or not 2 <= self.spatial_frames <= 201:
            raise ValueError("spatial_frames must be an integer between 2 and 201.")
        if self.event_log not in ("compact", "all"):
            raise ValueError("event_log must be compact or all.")


def signal_mobility(
    polymers: int | NDArray[np.int64], centre: float = 8.0,
) -> float | NDArray[np.float64]:
    """Evaluate S22 f(P) stably: (tanh((lambda-P)/2) + 1.2) / 2.2.

    tanh(z/2) equals (exp(z)-1)/(exp(z)+1), avoiding overflow in exp.
    P is the local polymer entity count, our stated lattice-concentration proxy.
    """
    return (np.tanh((centre - np.asarray(polymers)) / 2.0) + 1.2) / 2.2


def clusters(state: BindingState) -> tuple[tuple[int, ...], ...]:
    """Return connected bacterial components, including unbridged single cells."""
    parents = {b.identifier: b.identifier for b in state.bacteria}

    def find(identifier: int) -> int:
        while parents[identifier] != identifier:
            identifier = parents[identifier]
        return identifier

    # One array mask replaces one tiny NumPy reduction per polymer. Surviving
    # bridge rows retain their original order, including union-find tie choices.
    for contacts in state.polymer_contacts[np.all(state.polymer_contacts >= 0, axis=1)]:
        a, b = (int(state.occupancy.ravel()[site]) for site in contacts)
        parents[find(a)] = find(b)
    groups: dict[int, list[int]] = {}
    for identifier in parents:
        groups.setdefault(find(identifier), []).append(identifier)
    return tuple(sorted(tuple(sorted(members)) for members in groups.values()))


def cluster_can_move(state: BindingState, members: tuple[int, ...], direction: int) -> bool:
    """Allow a one-site translation only if every footprint remains unobstructed."""
    dy, dx = DIRECTIONS[direction]
    size = state.occupancy.shape[0]
    for bacterium in state.bacteria:
        if bacterium.identifier not in members:
            continue
        y, x = bacterium.y + dy, bacterium.x + dx
        if not (0 <= y < size - 1 and 0 <= x < size - 1):
            return False
        target = state.occupancy[y:y + 2, x:x + 2]
        if np.any((target >= 0) & ~np.isin(target, members)):
            return False
    return True


def retether(state: BindingState) -> None:
    """Anchor attached polymers at their first surviving contact (MOVE-03).

    A released last end leaves the polymer at its last contact. Bound signals
    remain with their polymer; no signal changes pool during reanchoring.
    """
    for polymer, contacts in enumerate(state.polymer_contacts):
        used = contacts[contacts >= 0]
        if len(used):
            state.polymer_sites[polymer] = used[0]


def move_cluster(
    state: BindingState, members: tuple[int, ...], direction: int, advection: bool = True,
) -> dict[str, int | str]:
    """Translate connected bacteria, their attachments, and local small entities.

    MOVE-04 shifts free chemicals in the old/new footprint union in the travel
    direction. At a closed outer edge a pushed chemical stays at the edge.
    Updates use the pre-move arrays, preventing double movement or count loss.
    """
    if not cluster_can_move(state, members, direction):
        raise ValueError("Cluster translation is blocked.")
    size = state.occupancy.shape[0]
    dy, dx = DIRECTIONS[direction]
    offset = dy * size + dx
    old = np.flatnonzero(np.isin(state.occupancy.ravel(), members))
    new = old + offset
    signal_count = polymer_count = 0
    unattached = np.all(state.polymer_contacts < 0, axis=1)
    # Determine attachment membership before changing the occupancy grid.
    attached = np.any(np.isin(state.polymer_contacts, old), axis=1)
    if advection:
        swept = np.union1d(old, new)
        y, x = np.divmod(swept, size)
        destinations = np.clip(y + dy, 0, size - 1) * size + np.clip(x + dx, 0, size - 1)
        counts = state.free_signals.ravel()[swept].copy()
        signal_count = int(counts.sum())
        state.free_signals.ravel()[swept] = 0
        np.add.at(state.free_signals.ravel(), destinations, counts)
        pushed = unattached & np.isin(state.polymer_sites, swept)
        polymer_count = int(np.count_nonzero(pushed))
        y, x = np.divmod(state.polymer_sites[pushed], size)
        state.polymer_sites[pushed] = (np.clip(y + dy, 0, size - 1) * size
                                       + np.clip(x + dx, 0, size - 1))
    # Move site-owned quantities simultaneously, including receptor-bound signals.
    for array, empty in ((state.occupancy, -1), (state.bacteria_signals, 0),
                         (state.bacteria_polymers, 0)):
        values = array.ravel()[old].copy()
        array.ravel()[old] = empty
        array.ravel()[new] = values
    contacts = state.polymer_contacts[attached].copy()
    contacts[contacts >= 0] += offset
    state.polymer_contacts[attached] = contacts
    state.polymer_sites[attached] += offset
    before = next(b for b in state.bacteria if b.identifier == members[0])
    state.bacteria = tuple(
        Bacterium(b.identifier, b.x + dx, b.y + dy) if b.identifier in members else b
        for b in state.bacteria
    )
    return {"reaction": "cluster_move", "x": before.x, "y": before.y,
            "to_x": before.x + dx, "to_y": before.y + dy, "polymer": -1,
            "bacterium": members[0], "cluster_size": len(members),
            "cluster_members": ",".join(map(str, members)),
            "advected_signals": signal_count, "advected_polymers": polymer_count}


@dataclass
class SpatialFrame:
    """One observation for maps and animation; arrays are copied when recorded."""

    time: float
    free_signals: NDArray[np.int64]
    bacteria: tuple[Bacterium, ...]
    sensed: NDArray[np.int64] | None = None
    expression: NDArray[np.float64] | None = None
    luminescence: NDArray[np.float64] | None = None


@dataclass
class MovementResult(BindingResult):
    """Add total event counts and bounded spatial frames to the common result."""

    total_events: int = 0
    event_counts: dict[str, int] = field(default_factory=dict)
    frames: list[SpatialFrame] = field(default_factory=list)
    positions: list[dict[str, int | float]] = field(default_factory=list)
    response_trajectory: list[dict[str, int | float]] = field(default_factory=list)
    cell_response: list[dict[str, int | float]] = field(default_factory=list)
    response_final: dict[str, NDArray] | None = None
    population_trajectory: list[dict[str, int | float]] = field(default_factory=list)
    lineage: list[dict[str, int | float]] = field(default_factory=list)
    produced_signals: int = 0
    scenario_note: str = ""


class MovementEngine:
    """Exact hierarchical SSA channels with local propensity updates.

    Signal hops dominate runtime. A hop changes only its source and destination;
    refreshing those sites avoids rebuilding thousands of unrelated channels.
    PB contacts refresh their dependent channels; actual translations/divisions
    rebuild geometry. Canonical partial sums retain the prior seeded arithmetic.
    The reference chemical solver remains available in binding mode.
    """

    def __init__(
        self, state: BindingState, binding: BindingConfig, movement: MovementConfig,
    ) -> None:
        self.state, self.binding, self.movement = state, binding, movement
        size = state.occupancy.shape[0]
        sites = np.arange(size * size)
        y, x = np.divmod(sites, size)
        self.destinations = np.full((size * size, 4), -1, dtype=np.int64)
        for direction, (dy, dx) in enumerate(DIRECTIONS):
            valid = (y + dy >= 0) & (y + dy < size) & (x + dx >= 0) & (x + dx < size)
            self.destinations[valid, direction] = sites[valid] + dy * size + dx
        self.rebuild()

    def rebuild(self) -> None:
        """Recompute all channels after a geometry/connectivity change."""
        state, config = self.state, self.movement
        nsites, npolymers = state.occupancy.size, len(state.polymer_sites)
        self.neighbours = neighbour_contacts(state)
        self.polymer_bins: list[set[int]] = [set() for _ in range(nsites)]
        for polymer, site in enumerate(state.polymer_sites):
            self.polymer_bins[int(site)].add(polymer)
        self.edge_factors = np.zeros((nsites, 4))
        occupied = state.occupancy.ravel()
        for direction in range(4):
            destination = self.destinations[:, direction]
            valid = destination >= 0
            obstructed = ((occupied[valid] >= 0) | (occupied[destination[valid]] >= 0))
            self.edge_factors[valid, direction] = np.where(
                obstructed, config.obstruction_factor, 1.0)
        self.bs_channels: list[list[Reaction]] = [[] for _ in range(nsites)]
        self.polymer_channels: list[list[Reaction]] = [[] for _ in range(npolymers)]
        self.signal_rates = np.zeros((nsites, 4))
        self.polymer_rates = np.zeros((npolymers, 4))
        self.bs_tree = RateTree(np.zeros(nsites))
        self.chem_tree = RateTree(np.zeros(npolymers))
        self.signal_tree = RateTree(np.zeros(nsites))
        self.polymer_tree = RateTree(np.zeros(npolymers))
        for site in range(nsites):
            self.refresh_site(site)
        for polymer in range(npolymers):
            self.refresh_polymer_motion(polymer)
        self.groups = clusters(state)
        self.cluster_rates = np.array([
            [config.d_b / len(members) if cluster_can_move(state, members, direction) else 0.0
             for direction in range(4)] for members in self.groups
        ], dtype=float).reshape(-1, 4)
        self.cluster_tree = RateTree(self.cluster_rates.sum(axis=1))
        # A binding target can be at capacity and have no positive PB channel.
        # Watch structural eligibility too, so releasing a full site restores
        # every newly available channel, including those absent from the cache.
        self.pb_watchers = [set() for _ in range(nsites)]
        self.pb_targets = [self._eligible_pb_targets(p) for p in range(npolymers)]
        for polymer, targets in enumerate(self.pb_targets):
            for site in targets:
                self.pb_watchers[site].add(polymer)
        # Stable sorting gives the old site-ascending, polymer-ID-ascending
        # insertion order. Recompute only when an anchor actually moves.
        self.chemical_order = np.argsort(state.polymer_sites, kind="stable")
        self.chemical_order_dirty = False

    def _eligible_pb_targets(self, polymer: int) -> tuple[int, ...]:
        """Index possible PB targets independently of their current propensity."""
        state = self.state
        first, second = map(int, state.polymer_contacts[polymer])
        if first < 0 and second < 0:
            site = int(state.polymer_sites[polymer])
            return (site,) if state.occupancy.ravel()[site] >= 0 else ()
        if first >= 0 and second >= 0:
            return ()
        contact = first if first >= 0 else second
        host = int(state.occupancy.ravel()[contact])
        return self.neighbours[host]

    def _refresh_polymer_chemistry(self, polymer: int) -> None:
        """Refresh one carrier after its signals or an eligible PB target changes."""
        local = len(self.polymer_bins[int(self.state.polymer_sites[polymer])])
        channels = polymer_reactions(self.state, self.binding, self.neighbours, polymer, local)
        self.polymer_channels[polymer] = channels
        self.chem_tree.update(polymer, math.fsum(c.propensity for c in channels))

    def _refresh_pb_index(self, polymer: int) -> None:
        """Maintain structural PB dependencies after contacts or anchors change."""
        for site in self.pb_targets[polymer]:
            self.pb_watchers[site].remove(polymer)
        self.pb_targets[polymer] = self._eligible_pb_targets(polymer)
        for site in self.pb_targets[polymer]:
            self.pb_watchers[site].add(polymer)

    def refresh_after_pb(self, polymer: int, target_site: int, kind: str) -> None:
        """Update one PB event without recomputing the entire physical geometry.

        The changed site's capacity affects all structurally eligible polymers.
        Reanchoring can also change two sites' PS accessibility and signal
        mobility. A second-end change may merge or split bacterial components.
        No event rates, channel order or random draws are approximated.
        """
        state = self.state
        affected = set(self.pb_watchers[target_site])
        affected.add(polymer)
        old_site = int(state.polymer_sites[polymer])
        first, second = map(int, state.polymer_contacts[polymer])
        new_site = first if first >= 0 else second if second >= 0 else old_site
        if new_site != old_site:
            self.polymer_bins[old_site].remove(polymer)
            self.polymer_bins[new_site].add(polymer)
            state.polymer_sites[polymer] = new_site
            self.chemical_order_dirty = True
        self._refresh_pb_index(polymer)

        refreshed = set()
        for site in sorted({old_site, new_site}):
            self.refresh_site(site)
            refreshed.update(self.polymer_bins[site])
        for other in sorted(affected - refreshed):
            self._refresh_polymer_chemistry(other)
        self.refresh_polymer_motion(polymer)

        ends = int(first >= 0) + int(second >= 0)
        bridge_changed = ends == 2 if kind == "pb_bind" else ends == 1
        if bridge_changed:
            groups = clusters(state)
            if groups != self.groups:
                self.groups = groups
                self.cluster_rates = np.array([
                    [self.movement.d_b / len(members)
                     if cluster_can_move(state, members, direction) else 0.0
                     for direction in range(4)] for members in groups
                ], dtype=float).reshape(-1, 4)
                self.cluster_tree = RateTree(self.cluster_rates.sum(axis=1))

        # v0.10.0 rebuilt zero trees after every PB event. Retaining incremental
        # sums would accumulate different floating-point rounding. Recreate the
        # exact old insertion order in NumPy, using cached unchanged leaf rates.
        self.bs_tree = RateTree.from_sequential_updates(self.bs_tree.values)
        self.signal_tree = RateTree.from_sequential_updates(self.signal_tree.values)
        if self.chemical_order_dirty:
            self.chemical_order = np.argsort(state.polymer_sites, kind="stable")
            self.chemical_order_dirty = False
        self.chem_tree = RateTree.from_sequential_updates(self.chem_tree.values,
                                                         self.chemical_order)
        self.polymer_tree = RateTree.from_sequential_updates(self.polymer_tree.values)

    @property
    def trees(self) -> tuple[RateTree, ...]:
        return (self.bs_tree, self.chem_tree, self.signal_tree,
                self.polymer_tree, self.cluster_tree)

    def reset_partial_sums(self) -> None:
        """Rebuild weighted sums from cached rates without recalculating geometry."""
        self.bs_tree = RateTree(self.bs_tree.values)
        self.chem_tree = RateTree(self.chem_tree.values)
        self.signal_tree = RateTree(self.signal_tree.values)
        self.polymer_tree = RateTree(self.polymer_tree.values)
        self.cluster_tree = RateTree(self.cluster_tree.values)

    def refresh_site(self, site: int) -> None:
        """Refresh channels affected by a local free-signal or polymer count change."""
        state = self.state
        signals = int(state.free_signals.ravel()[site])
        bound = int(state.bacteria_signals.ravel()[site])
        if signals < 0 or not 0 <= bound <= BS_PER_SITE:
            raise ValueError("A local signal count or receptor capacity was violated.")
        channels = bacterial_reactions(state, self.binding, site)
        self.bs_channels[site] = channels
        self.bs_tree.update(site, math.fsum(c.propensity for c in channels))
        local = len(self.polymer_bins[site])
        mobility = float(signal_mobility(local, self.movement.crowding_lambda))
        self.signal_rates[site] = signals * self.movement.d_s * mobility * self.edge_factors[site]
        self.signal_tree.update(site, float(self.signal_rates[site].sum()))
        for polymer in sorted(self.polymer_bins[site]):
            channels = polymer_reactions(state, self.binding, self.neighbours, polymer, local)
            self.polymer_channels[polymer] = channels
            self.chem_tree.update(polymer, math.fsum(c.propensity for c in channels))

    def refresh_polymer_motion(self, polymer: int) -> None:
        """Only polymers without bacterial contacts diffuse independently."""
        site = int(self.state.polymer_sites[polymer])
        free = np.all(self.state.polymer_contacts[polymer] < 0)
        self.polymer_rates[polymer] = (self.movement.d_p * self.edge_factors[site]
                                        if free else 0.0)
        self.polymer_tree.update(polymer, float(self.polymer_rates[polymer].sum()))

    @staticmethod
    def choose(weights: list[float] | NDArray[np.float64], rng: np.random.Generator) -> int:
        """Select a small within-family channel by its conditional probability."""
        cumulative = np.cumsum(weights)
        draw = float(rng.random()) * float(cumulative[-1])
        return min(int(np.searchsorted(cumulative, draw, side="right")), len(weights) - 1)

    def event(self, family: int, index: int, rng: np.random.Generator) -> dict[str, int | str]:
        """Apply a selected channel and update exactly the affected rate caches."""
        state = self.state
        size = state.occupancy.shape[0]
        if family in (0, 1):
            channels = self.bs_channels[index] if family == 0 else self.polymer_channels[index]
            reaction = channels[self.choose([c.propensity for c in channels], rng)]
            record = apply_reaction(state, reaction)
            if reaction.kind.startswith("pb_"):
                site = int(record["y"]) * size + int(record["x"])
                self.refresh_after_pb(reaction.target, site, reaction.kind)
            else:
                self.refresh_site(int(record["y"]) * size + int(record["x"]))
            return record
        if family == 4:
            direction = self.choose(self.cluster_rates[index], rng)
            record = move_cluster(state, self.groups[index], direction, self.movement.advection)
            self.rebuild()
            return record
        rates = self.signal_rates[index] if family == 2 else self.polymer_rates[index]
        direction = self.choose(rates, rng)
        source = index if family == 2 else int(state.polymer_sites[index])
        destination = int(self.destinations[source, direction])
        if destination < 0:
            raise ValueError("A movement channel crossed a closed boundary.")
        if family == 2:
            state.free_signals.ravel()[source] -= 1
            state.free_signals.ravel()[destination] += 1
        else:
            self.polymer_bins[source].remove(index)
            self.polymer_bins[destination].add(index)
            state.polymer_sites[index] = destination
            self._refresh_pb_index(index)
            self.chemical_order_dirty = True
            self.refresh_polymer_motion(index)
        self.refresh_site(source)
        self.refresh_site(destination)
        sy, sx = divmod(source, size)
        ty, tx = divmod(destination, size)
        return {"reaction": "signal_hop" if family == 2 else "polymer_hop",
                "x": sx, "y": sy, "to_x": tx, "to_y": ty,
                "polymer": -1 if family == 2 else index, "bacterium": -1}


def simulate_movement(
    initial: InitialState, binding: BindingConfig, movement: MovementConfig, seed: int,
    response: ResponseConfig | None = None, strain: str = "MM32",
    population: PopulationConfig | None = None,
    progress: Callable[[float, int], None] | None = None,
) -> MovementResult:
    """Sample coupled chemistry/transport without fixed time steps or tau leaping.

    Optional progress observes actual model time every 1000 events and at the
    endpoints. It does not schedule events, sample the model, or draw randomness.
    """
    if type(seed) is not int or not 0 <= seed < 2**64:
        raise ValueError("seed must be an integer between 0 and 2**64 - 1.")
    total_signals = int(initial.free_signals.sum())
    total_polymers = int(initial.free_polymers.sum())
    binding.validate_polymer_rates(total_polymers)
    state = BindingState.from_initial(initial)
    state.validate(total_signals, total_polymers)
    tracker = ResponseTracker(state, response, strain) if response is not None else None
    engine = MovementEngine(state, binding, movement)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([seed, 2])))
    result = MovementResult(state, "complete", 0.0)
    # Independent streams keep lifetime and placement choices separate from SSA
    # events. Omitting population preserves the earlier diagnostic event history.
    schedule = None
    production_rate = 0.0
    placement_rng = None
    if population is not None:
        lifetime_rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([seed, 3])))
        placement_rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([seed, 4])))
        schedule = DivisionSchedule(state, population, lifetime_rng)
        production_rate = population.signal_rate(strain)
        result.lineage = schedule.lineage
        if population.growth and (population.doubling_time != 15000.0
                                  or population.division_width != 20000.0):
            result.scenario_note = "Division timing differs from source"
    counts: Counter[str] = Counter()
    sample_index = frame_index = 0
    frame_times = np.linspace(0.0, binding.duration, movement.spatial_frames)

    def record_sample(time: float) -> None:
        state.validate(total_signals, total_polymers)
        result.trajectory.append(sample_counts(state, time))
        if population is not None:
            result.population_trajectory.append({
                "time": time, "bacteria": len(state.bacteria),
                "divisions": counts["division"], "blocked_attempts": counts["division_blocked"],
                "produced_signals": result.produced_signals,
                "expected_total_signals": total_signals,
                "actual_total_signals": int(state.free_signals.sum()
                                            + state.bacteria_signals.sum()
                                            + state.polymer_signals.sum()),
            })
        result.positions.extend({"time": time, "bacterium": b.identifier, "x": b.x, "y": b.y}
                                for b in state.bacteria)
        if tracker is not None:
            tracker.validate(state, time)
            values = tracker.observe(time)
            result.response_trajectory.append(response_summary(values, time))
            for index, bacterium in enumerate(state.bacteria):
                result.cell_response.append({
                    "time": time, "bacterium": bacterium.identifier,
                    "x": bacterium.x, "y": bacterium.y,
                    "sensed": int(values["sensed"][index]),
                    "expression": float(values["expression"][index]),
                    "luminescence": float(values["luminescence"][index]),
                })

    def record_frame(time: float) -> None:
        values = tracker.observe(time) if tracker is not None else None
        result.frames.append(SpatialFrame(
            time, state.free_signals.copy(), state.bacteria,
            values["sensed"] if values is not None else None,
            values["expression"] if values is not None else None,
            values["luminescence"] if values is not None else None,
        ))

    def observations(limit: float, inclusive: bool) -> None:
        nonlocal sample_index, frame_index
        def due(time: float) -> bool:
            return time <= limit if inclusive else time < limit
        while sample_index * binding.sample_interval <= binding.duration:
            time = sample_index * binding.sample_interval
            if not due(time):
                break
            record_sample(time)
            sample_index += 1
        while frame_index < len(frame_times) and due(float(frame_times[frame_index])):
            record_frame(float(frame_times[frame_index]))
            frame_index += 1

    observations(0, True)
    if progress is not None:
        progress(0.0, 0)
    while True:
        totals = [tree.total for tree in engine.trees]
        if population is not None:
            # Equal per-cell production rates give the total population rate.
            # The producer and its footprint site are selected only if it fires.
            totals.append(production_rate * len(state.bacteria))
        total = math.fsum(totals)
        if not math.isfinite(total) or total < 0:
            raise ValueError("Total event rate is invalid; check configuration magnitudes.")
        stochastic_time = result.time + float(rng.exponential(1.0 / total)) if total else math.inf
        division_time = schedule.next_time if schedule is not None else math.inf
        is_division = division_time <= stochastic_time and math.isfinite(division_time)
        next_time = division_time if is_division else stochastic_time
        if next_time > binding.duration:
            observations(binding.duration, True)
            result.time = binding.duration
            break
        if next_time < result.time or (next_time == result.time and not is_division):
            raise ArithmeticError("Event time cannot advance at this numerical precision.")
        observations(next_time, False)
        if is_division:
            # Exponential waiting times may be resampled after this appointment:
            # between discrete changes they are memoryless. Triangular lifetimes
            # themselves are never replaced by a constant division hazard.
            record = schedule.divide(state, next_time, placement_rng)
            if record["reaction"] == "division":
                if tracker is not None:
                    tracker.add_daughter(int(record["daughter"]), next_time,
                                         population.daughter_expression)
                engine.rebuild()  # new footprint, receptors, obstruction, and singleton cluster
        else:
            cumulative = np.cumsum(totals)
            draw = float(rng.random()) * float(cumulative[-1])
            family = min(int(np.searchsorted(cumulative, draw, side="right")), len(totals) - 1)
            if family == len(engine.trees):
                record = produce_signal(state, rng)
                result.produced_signals += 1
                total_signals += 1  # production changes the conservation target by exactly one
                engine.refresh_site(int(record["y"]) * state.occupancy.shape[0] + int(record["x"]))
            else:
                weight = draw - (float(cumulative[family - 1]) if family else 0.0)
                tree = engine.trees[family]
                # Guard a final-ULP rounding to total without selecting a zero channel.
                weight = max(0.0, min(weight, float(np.nextafter(tree.total, 0.0))))
                index = tree.select(weight)
                record = engine.event(family, index, rng)
        if tracker is not None and record["reaction"] in ("bs_bind", "bs_unbind"):
            tracker.change_sensing(int(record["bacterium"]),
                                   1 if record["reaction"] == "bs_bind" else -1, next_time)
        result.time = next_time
        result.total_events += 1
        if progress is not None and result.total_events % 1000 == 0:
            progress(result.time, result.total_events)
        kind = str(record["reaction"])
        counts[kind] += 1
        if movement.event_log == "all" or kind != "signal_hop":
            result.events.append({"event": result.total_events, "time": next_time, **record})
        if result.total_events % movement.validate_every == 0:
            state.validate(total_signals, total_polymers)
            if tracker is not None:
                tracker.validate(state, result.time)
            # Periodic reconstruction also prevents long-run partial-sum drift.
            engine.reset_partial_sums()
        # With zero-width lifetimes several appointments can coincide. Sample
        # the state after all appointments at that timestamp (unless event-capped).
        simultaneous = schedule is not None and schedule.next_time == next_time
        if not simultaneous:
            observations(next_time, True)
        if result.total_events >= binding.max_events and (
            result.time < binding.duration or simultaneous
        ):
            result.status = "event_limit"
            break
    if result.trajectory[-1]["time"] != result.time:
        record_sample(result.time)
    if result.frames[-1].time != result.time:
        record_frame(result.time)
    result.event_counts = dict(counts)
    state.validate(total_signals, total_polymers)
    if tracker is not None:
        tracker.validate(state, result.time)
        result.response_final = tracker.observe(result.time)
    if progress is not None:
        progress(result.time, result.total_events)
    return result
