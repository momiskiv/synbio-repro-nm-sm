# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Transport geometry, carrier conservation, cache correctness, and diffusion law."""

import copy
import unittest

import numpy as np

from polymer_clustering.binding import BindingConfig, BindingState, apply_reaction, reaction_channels
from polymer_clustering.model import Bacterium, InitialConfig, InitialState, initialise
from polymer_clustering.movement import (
    MovementConfig, MovementEngine, cluster_can_move, clusters, move_cluster,
    retether, signal_mobility, simulate_movement,
)
from polymer_clustering.rates import RateTree
from test_binding import channel, rates, small_world


def empty_world(size: int = 5, signals: int = 1) -> InitialState:
    """Put signals at the centre of an otherwise empty closed square lattice."""
    free = np.zeros((size, size), dtype=np.int64)
    free[size // 2, size // 2] = signals
    return InitialState(np.full((size, size), -1, dtype=np.int32), (), free, np.zeros_like(free))


class RateTreeTests(unittest.TestCase):
    def test_weighted_intervals_zeros_and_updates(self) -> None:
        tree = RateTree(np.array([0.0, 2.0, 0.0, 3.0]))
        self.assertEqual(tree.total, 5)
        self.assertEqual([tree.select(x) for x in (0, 1.9, 2, 4.9)], [1, 1, 3, 3])
        tree.update(1, 0)
        self.assertEqual(tree.select(0), 3)
        tree.update(3, 0)
        self.assertEqual(tree.total, 0)
        with self.assertRaises(ValueError):
            tree.select(0)
        with self.assertRaises(ValueError):
            tree.update(0, -1)

    def test_updates_match_independent_cumulative_sums(self) -> None:
        rng = np.random.default_rng(4)
        values = rng.random(127)
        tree = RateTree(values)
        for _ in range(500):
            index, value = int(rng.integers(len(values))), float(rng.random())
            values[index] = value
            tree.update(index, value)
            self.assertAlmostEqual(tree.total, float(values.sum()), places=10)
            draw = float(rng.random()) * float(values.sum())
            self.assertEqual(tree.select(draw), int(np.searchsorted(np.cumsum(values), draw, side="right")))


class MovementTests(unittest.TestCase):
    def test_independent_placement_and_legacy_replay(self) -> None:
        without = initialise(InitialConfig(polymers=0))
        with_polymer = initialise(InitialConfig(polymers=100))
        self.assertEqual(without.bacteria, with_polymer.bacteria)
        np.testing.assert_array_equal(without.free_signals, with_polymer.free_signals)
        # This v0.2 seed had 382 signals initially sharing bacterial sites.
        legacy = initialise(InitialConfig(rng_scheme="legacy_v02"))
        self.assertEqual(int(legacy.free_signals[legacy.occupancy >= 0].sum()), 382)
        self.assertFalse(np.array_equal(legacy.free_signals, without.free_signals))

    def test_stable_published_mobility_function(self) -> None:
        values = signal_mobility(np.array([0, 8, 100000], dtype=np.int64))
        self.assertTrue(np.all(np.isfinite(values)))
        self.assertTrue(np.all(np.diff(values) < 0))
        self.assertAlmostEqual(values[0], ((np.exp(8) - 1) / (np.exp(8) + 1) + 1.2) / 2.2)
        self.assertAlmostEqual(values[1], 1.2 / 2.2)
        self.assertAlmostEqual(values[2], 0.2 / 2.2)

    def test_directional_signal_rates_boundaries_and_obstruction(self) -> None:
        state = BindingState.from_initial(small_world(polymers=0, signals=4))
        engine = MovementEngine(state, BindingConfig(), MovementConfig())
        # At the upper-left array corner only +y and +x exist. These edges touch
        # a bacterial footprint and receive the explicitly chosen factor 0.1.
        expected = 4 * 0.05 * float(signal_mobility(0)) * 0.1
        np.testing.assert_allclose(engine.signal_rates[0], [0, expected, 0, expected])
        blocked = MovementEngine(state, BindingConfig(), MovementConfig(obstruction_factor=0))
        self.assertEqual(blocked.signal_rates[0].sum(), 0)
        open_engine = MovementEngine(state, BindingConfig(), MovementConfig(obstruction_factor=1))
        self.assertAlmostEqual(open_engine.signal_rates[0, 1], expected * 10)

    def test_signal_reaches_a_different_site_and_stays_inside(self) -> None:
        initial = empty_world(size=3, signals=10)
        result = simulate_movement(initial, BindingConfig(duration=10), MovementConfig(), 2)
        self.assertGreater(result.event_counts.get("signal_hop", 0), 0)
        self.assertEqual(result.state.free_signals.sum(), 10)
        self.assertLess(result.state.free_signals[1, 1], 10)
        self.assertEqual(initial.free_signals[1, 1], 10)

    def test_polymer_motion_carries_bound_signal(self) -> None:
        state = BindingState.from_initial(small_world())
        apply_reaction(state, channel(state, "ps_bind"))
        engine = MovementEngine(state, rates(), MovementConfig(d_s=0, d_b=0))
        engine.event(3, 0, np.random.default_rng(2))
        self.assertNotEqual(state.polymer_sites[0], 0)
        self.assertEqual(state.polymer_signals[0], 1)
        state.validate(20, 1)
        apply_reaction(state, channel(state, "pb_bind"))
        retether(state)
        engine.rebuild()
        self.assertEqual(engine.polymer_rates[0].sum(), 0)

    def test_bridge_moves_as_one_unit_and_retains_attachments(self) -> None:
        state = BindingState.from_initial(small_world())
        apply_reaction(state, channel(state, "pb_bind"))
        apply_reaction(state, channel(state, "pb_bind"))
        apply_reaction(state, channel(state, "ps_bind"))
        apply_reaction(state, channel(state, "bs_bind"))
        self.assertEqual(clusters(state), ((0, 1),))
        before_contacts = state.polymer_contacts.copy()
        record = move_cluster(state, (0, 1), 1)
        state.validate(20, 1)
        self.assertEqual(record["cluster_size"], 2)
        self.assertEqual([(b.x, b.y) for b in state.bacteria], [(0, 1), (2, 1)])
        np.testing.assert_array_equal(state.polymer_contacts, before_contacts + 4)
        self.assertEqual(state.polymer_sites[0], 4)
        self.assertEqual(state.polymer_signals[0], 1)
        self.assertEqual(state.bacteria_signals[1, 0], 1)
        self.assertEqual(clusters(state), ((0, 1),))

    def test_cluster_inverse_size_rate_and_blocking(self) -> None:
        state = BindingState.from_initial(small_world())
        self.assertFalse(cluster_can_move(state, (0,), 3))  # Another cell is to its right.
        self.assertFalse(cluster_can_move(state, (0,), 0))  # Closed lower edge.
        apply_reaction(state, channel(state, "pb_bind"))
        apply_reaction(state, channel(state, "pb_bind"))
        engine = MovementEngine(state, rates(), MovementConfig())
        self.assertAlmostEqual(engine.cluster_rates[0, 1], 0.0005 / 2)
        self.assertEqual(engine.cluster_rates[0, 0], 0)
        self.assertEqual(engine.cluster_rates[0, 2], 0)
        self.assertEqual(engine.cluster_rates[0, 3], 0)

    def test_bridge_release_splits_components_and_reanchors_surviving_end(self) -> None:
        state = BindingState.from_initial(small_world())
        apply_reaction(state, channel(state, "pb_bind"))
        apply_reaction(state, channel(state, "pb_bind"))
        apply_reaction(state, channel(state, "pb_unbind"))
        retether(state)
        self.assertEqual(clusters(state), ((0,), (1,)))
        self.assertEqual(state.polymer_sites[0], state.polymer_contacts[0, 1])
        state.validate(20, 1)

    def test_advection_conserves_and_displaces_front_chemicals(self) -> None:
        occupancy = np.full((5, 5), -1, dtype=np.int32)
        occupancy[1:3, 1:3] = 0
        free = np.zeros((5, 5), dtype=np.int64)
        free[1, 1], free[1, 3] = 3, 4
        polymers = np.zeros_like(free)
        polymers[1, 3] = 1
        initial = InitialState(occupancy, (Bacterium(0, 1, 1),), free, polymers)
        state = BindingState.from_initial(initial)
        record = move_cluster(state, (0,), 3)
        state.validate(7, 1)
        self.assertEqual(state.free_signals[1, 2], 3)
        self.assertEqual(state.free_signals[1, 4], 4)
        self.assertEqual(state.polymer_sites[0], 9)
        self.assertEqual(record["advected_signals"], 7)
        self.assertEqual(record["advected_polymers"], 1)
        # Pushing into an outer boundary retains the chemical inside the grid.
        move_cluster(state, (0,), 3)
        state.validate(7, 1)
        self.assertEqual(state.free_signals[1, 4], 4)
        self.assertEqual(state.polymer_sites[0], 9)

    def test_no_advection_leaves_free_chemicals_at_their_sites(self) -> None:
        state = BindingState.from_initial(small_world(polymers=0))
        before = state.free_signals.copy()
        move_cluster(state, (0,), 1, advection=False)
        np.testing.assert_array_equal(state.free_signals, before)
        state.validate(20, 0)

    def test_cached_channels_match_full_reconstruction_after_events(self) -> None:
        state = BindingState.from_initial(small_world(polymers=3, signals=30))
        config = rates(k_pb_on=0.1, k_pb_off=0.1, k_ps_on=0.1, k_ps_off=0.1)
        movement = MovementConfig(d_s=0.3, d_p=0.5, d_b=0.2)
        engine = MovementEngine(state, config, movement)
        rng = np.random.default_rng(7)
        # Force a mobile polymer event before strong bacterial capture removes
        # that opportunity; random coverage alone is not a reliable assertion.
        kinds = {engine.event(3, 0, rng)["reaction"]}
        state.validate(30, 3)
        reference = MovementEngine(copy.deepcopy(state), config, movement)
        for actual, expected in zip(engine.trees, reference.trees):
            np.testing.assert_allclose(actual.values, expected.values, rtol=1e-12, atol=1e-12)
        for _ in range(200):
            rates_by_family = [tree.total for tree in engine.trees]
            family = engine.choose(rates_by_family, rng)
            tree = engine.trees[family]
            index = tree.select(float(rng.random()) * tree.total)
            kinds.add(engine.event(family, index, rng)["reaction"])
            state.validate(30, 3)
            reference = MovementEngine(copy.deepcopy(state), config, movement)
            for actual, expected in zip(engine.trees, reference.trees):
                np.testing.assert_allclose(actual.values, expected.values, rtol=1e-12, atol=1e-12)
            actual_chem = [c for bucket in engine.bs_channels + engine.polymer_channels for c in bucket]
            expected_chem = reaction_channels(state, config, reference.neighbours)
            key = lambda c: (c.kind, c.target, c.contact, c.propensity)
            self.assertEqual(sorted(map(key, actual_chem)), sorted(map(key, expected_chem)))
        self.assertIn("cluster_move", kinds)
        self.assertIn("polymer_hop", kinds)
        self.assertIn("signal_hop", kinds)

    def test_sampling_logging_and_frame_controls_do_not_change_final_state(self) -> None:
        initial = small_world(polymers=2)
        first = simulate_movement(initial, rates(duration=5, sample_interval=1),
                                  MovementConfig(event_log="all", spatial_frames=3), 6)
        second = simulate_movement(initial, rates(duration=5, sample_interval=0.3),
                                   MovementConfig(event_log="compact", spatial_frames=8), 6)
        self.assertEqual(first.total_events, second.total_events)
        self.assertEqual(first.event_counts, second.event_counts)
        np.testing.assert_array_equal(first.state.free_signals, second.state.free_signals)
        np.testing.assert_array_equal(first.state.polymer_sites, second.state.polymer_sites)
        self.assertEqual(first.state.bacteria, second.state.bacteria)
        self.assertEqual(first.frames[-1].time, 5)
        self.assertEqual(second.frames[-1].time, 5)

    def test_event_cap_and_zero_rates_preserve_real_stop_time(self) -> None:
        cap = simulate_movement(empty_world(), BindingConfig(duration=10, max_events=1),
                                MovementConfig(), 1)
        self.assertEqual(cap.status, "event_limit")
        self.assertEqual(cap.frames[-1].time, cap.time)
        self.assertEqual(cap.trajectory[-1]["time"], cap.time)
        self.assertLess(cap.time, 10)
        stopped = simulate_movement(empty_world(), BindingConfig(duration=10),
                                    MovementConfig(d_s=0, d_p=0, d_b=0), 1)
        self.assertEqual(stopped.total_events, 0)
        self.assertEqual(stopped.time, 10)

    def test_free_diffusion_mean_squared_displacement(self) -> None:
        # In a large unobstructed 2D lattice, each direction has rate D*f(0).
        # One particle therefore has E[hops] = E[r^2] = 4*D*f(0)*t before edges.
        # 600 fixed seeds test the transport law independently of cache internals.
        hops, distances = [], []
        initial = empty_world(size=11)
        config = BindingConfig(duration=10, sample_interval=10)
        for seed in range(600):
            result = simulate_movement(initial, config, MovementConfig(spatial_frames=2), seed)
            y, x = np.argwhere(result.state.free_signals)[0]
            hops.append(result.event_counts.get("signal_hop", 0))
            distances.append(int((y - 5) ** 2 + (x - 5) ** 2))
        expected = 4 * 0.05 * float(signal_mobility(0)) * 10
        self.assertAlmostEqual(float(np.mean(hops)), expected, delta=0.2)
        self.assertAlmostEqual(float(np.mean(distances)), expected, delta=0.35)

    def test_invalid_movement_settings(self) -> None:
        for settings in ({"d_s": -1}, {"d_b": float("nan")}, {"obstruction_factor": 2},
                         {"advection": "yes"}, {"validate_every": 0},
                         {"spatial_frames": 202}, {"event_log": "unknown"}):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                MovementConfig(**settings)


if __name__ == "__main__":
    unittest.main()
