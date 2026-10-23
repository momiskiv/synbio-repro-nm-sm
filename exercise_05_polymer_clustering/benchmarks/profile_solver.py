# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Profile bounded software diagnostics without plotting or changing model rules.

Run this file from any directory. --source-root selects the extracted release
to measure; comparing its result signatures checks exact seeded replay. Timed
work excludes imports, rendering, profiling overhead and signature generation.
The forced PB cycle is a valid-channel cache benchmark, not a biological run.
"""

import argparse
import cProfile
from dataclasses import fields, is_dataclass
import hashlib
import json
from pathlib import Path
import pstats
import sys
import time

import numpy as np


def signature(value):
    """Hash every result field and array byte, including event times and identity."""
    if isinstance(value, np.ndarray):
        return dict(dtype=str(value.dtype), shape=value.shape,
                    sha256=hashlib.sha256(value.tobytes()).hexdigest())
    if is_dataclass(value):
        return {f.name: signature(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, dict):
        return {k: signature(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [signature(v) for v in value]
    return value


def digest(value):
    return hashlib.sha256(json.dumps(signature(value), sort_keys=True,
                                     allow_nan=False).encode()).hexdigest()


def crowded_binding_state(count):
    """Construct a valid late-density stress state; this is not a simulated result.

    Seventy-nine cells form four separated packs of at most twenty cells.
    Contacts respect 30 ends/site, bridges join adjacent different bacteria,
    signals total 10000 and each polymer retains at most five signals. The
    fixture exercises occupied PS slots, many bridges and local crowding.
    """
    from polymer_clustering.binding import BindingState
    from polymer_clustering.model import Bacterium, InitialConfig, initialise
    state = BindingState.from_initial(initialise(InitialConfig(
        bacteria=79, polymers=count, seed=2026, rng_scheme="independent_v1")))
    bacteria = []
    state.occupancy[:] = -1
    for oy, ox in ((2, 2), (2, 18), (18, 2), (18, 18)):
        for row in range(4):
            for column in range(5):
                if len(bacteria) == 79:
                    break
                y, x = oy + 2*row, ox + 2*column
                bacterium = Bacterium(len(bacteria), x, y)
                bacteria.append(bacterium)
                state.occupancy[y:y+2, x:x+2] = bacterium.identifier
    state.bacteria = tuple(bacteria)
    occupied = state.occupancy.ravel()
    edges = []
    for site in np.flatnonzero(occupied >= 0):
        y, x = divmod(int(site), 50)
        for target in ((int(site)+1 if x < 49 else -1),
                       (int(site)+50 if y < 49 else -1)):
            if target >= 0 and occupied[target] >= 0 and occupied[site] != occupied[target]:
                edges.append((int(site), target))
    capacity = state.bacteria_polymers.ravel()
    bridges = min(int(.44*count), 3000)
    attached = min(int(.61*count), 6000)
    pointer = 0
    for polymer in range(bridges):
        for _ in range(len(edges)):
            a, b = edges[pointer % len(edges)]
            pointer += 1
            if capacity[a] < 30 and capacity[b] < 30:
                break
        else:
            raise ValueError("Crowded fixture exhausted bridge capacity.")
        state.polymer_contacts[polymer] = (a, b)
        state.polymer_sites[polymer] = a
        capacity[a] += 1
        capacity[b] += 1
    sites = np.flatnonzero(occupied >= 0)
    pointer = 0
    for polymer in range(bridges, attached):
        while capacity[sites[pointer % len(sites)]] >= 30:
            pointer += 1
        site = int(sites[pointer % len(sites)])
        pointer += 1
        state.polymer_contacts[polymer, 0] = site
        state.polymer_sites[polymer] = site
        capacity[site] += 1
    bound = min(6300, 5*count)
    state.polymer_signals[:] = bound//count
    state.polymer_signals[:bound % count] += 1
    receptors = min(2370, 10000-bound)
    state.bacteria_signals.ravel()[sites] = receptors//len(sites)
    state.bacteria_signals.ravel()[sites[:receptors % len(sites)]] += 1
    free = 10000-bound-receptors
    state.free_signals[:] = free//state.occupancy.size
    state.free_signals.ravel()[:free % state.occupancy.size] += 1
    state.validate(10000, count)
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--polymers", type=int, nargs="+", default=[1000, 3000, 7000])
    parser.add_argument("--case", choices=("pb-cycle", "crowded-pb", "trajectory"), default="pb-cycle")
    parser.add_argument("--cycles", type=int, default=10)
    parser.add_argument("--events", type=int, default=3000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--output", type=Path, default=Path("results/exercise5/performance"))
    args = parser.parse_args()
    if min(args.cycles, args.events, args.repeats) < 1 or min(args.polymers) < 0:
        parser.error("Counts must be positive; polymer count may be zero for trajectories.")
    if args.case != "trajectory" and min(args.polymers) == 0:
        parser.error("PB cycles require polymers.")
    sys.path.insert(0, str(args.source_root.expanduser().resolve()))
    from polymer_clustering.binding import BindingConfig, BindingState
    from polymer_clustering.model import InitialConfig, initialise
    from polymer_clustering.movement import MovementConfig, MovementEngine, simulate_movement
    from polymer_clustering.output import create_run_directory, software_metadata, write_json
    from polymer_clustering.population import PopulationConfig
    from polymer_clustering.response import ResponseConfig

    output = create_run_directory(args.output, "full")
    report = dict(purpose="Bounded software diagnostics; not assignment-scale predictions.",
                  software=software_metadata(), case=args.case, repeats=args.repeats,
                  source_root=str(args.source_root.resolve()),
                  timing_scope="Unprofiled solver/cache work only; no plotting or file output.",
                  results=[])
    for count in args.polymers:
        initial = initialise(InitialConfig(polymers=count, seed=2026, rng_scheme="independent_v1"))
        binding = BindingConfig(duration=30000, sample_interval=100, max_events=args.events,
                                k_pb_on=.005, k_pb_off=.001, k_ps_on=.01, k_ps_off=.001)
        movement = MovementConfig(spatial_frames=2, event_log="all")

        def prepare():
            if args.case == "trajectory":
                return None
            state = (crowded_binding_state(count) if args.case == "crowded-pb"
                     else BindingState.from_initial(initial))
            engine = MovementEngine(state, binding, movement)
            attached_ends = np.count_nonzero(state.polymer_contacts >= 0, axis=1)
            candidates = np.flatnonzero(attached_ends == (1 if args.case == "crowded-pb" else 0))
            candidates = [int(p) for p in candidates
                          if any(c.kind == "pb_bind" for c in engine.polymer_channels[p])]
            if not len(candidates):
                raise ValueError("This seeded state has no polymer on a bacterial footprint.")
            return engine, int(candidates[0])

        def run(prepared):
            if args.case == "trajectory":
                return simulate_movement(initial, binding, movement, 2026,
                                         response=ResponseConfig(), population=PopulationConfig())
            engine, polymer = prepared
            # Select a valid PB channel by an interior conditional draw. Every
            # cycle binds and releases one end, restoring the same physical state.
            class Draw:
                def random(self):
                    return self.value
            draw = Draw()
            for _ in range(args.cycles):
                new_end = int(np.flatnonzero(engine.state.polymer_contacts[polymer] < 0)[0])
                for kind in ("pb_bind", "pb_unbind"):
                    channels = engine.polymer_channels[polymer]
                    index = next(i for i, c in enumerate(channels) if c.kind == kind
                                 and (kind == "pb_bind" or c.contact == new_end))
                    cumulative = np.cumsum([c.propensity for c in channels])
                    before = cumulative[index - 1] if index else 0.
                    draw.value = float((before + cumulative[index]) / 2 / cumulative[-1])
                    record = engine.event(1, polymer, draw)
                    assert record["reaction"] == kind
            engine.state.validate(10000, count)
            return dict(state=engine.state,
                        rate_values=[tree.values.copy() for tree in engine.trees],
                        partial_sums=[tree.tree.copy() for tree in engine.trees],
                        groups=engine.groups)

        wall, cpu, hashes = [], [], []
        for _ in range(args.repeats):
            prepared = prepare()
            wall_start, cpu_start = time.perf_counter(), time.process_time()
            result = run(prepared)
            wall.append(time.perf_counter() - wall_start)
            cpu.append(time.process_time() - cpu_start)
            hashes.append(digest(result))
        if len(set(hashes)) != 1:
            raise ValueError("Repeated identical diagnostics produced different results.")
        record = dict(polymers=count, seed=2026, wall_seconds=wall, cpu_seconds=cpu,
                      median_wall_seconds=float(np.median(wall)),
                      median_cpu_seconds=float(np.median(cpu)), result_sha256=hashes[0])
        if args.case == "trajectory":
            record.update(events=result.total_events, stop_time=result.time,
                          status=result.status, event_counts=result.event_counts)
        else:
            record.update(events=2*args.cycles, cycles=args.cycles)
        if args.profile:
            prepared = prepare()
            profiler = cProfile.Profile()
            profiled = profiler.runcall(run, prepared)
            assert digest(profiled) == hashes[0]
            profile_path = output/f"{args.case}_{count}.prof"
            profiler.dump_stats(profile_path)
            stats = pstats.Stats(profiler)
            top = sorted(stats.stats.items(), key=lambda x: x[1][3], reverse=True)[:20]
            record["profile_total_seconds"] = stats.total_tt
            record["profile_top_cumulative"] = [dict(
                file=Path(key[0]).name, line=key[1], function=key[2],
                calls=value[1], own_seconds=value[2], cumulative_seconds=value[3])
                for key, value in top]
        report["results"].append(record)
        write_json(output/"benchmark.json", report)
        print(f"{args.case}: P={count}; median CPU={record['median_cpu_seconds']:.6f}s; "
              f"signature={hashes[0]}", flush=True)
    print(f"Diagnostic outputs: {output}", flush=True)


if __name__ == "__main__":
    main()
