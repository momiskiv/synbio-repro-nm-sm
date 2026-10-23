# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Compare local PB updates with full reconstruction at difficult boundaries."""

import copy
import unittest
from unittest.mock import patch

import numpy as np

from polymer_clustering.binding import BindingState
from polymer_clustering.batch import (
    OPTIMISED_MOVEMENT_SHA256, OPTIMISED_RATES_SHA256, PRE_OPTIMISATION_RATES_SHA256,
    PROGRESS_MOVEMENT_SHA256, numerical_sources_match,
)
from polymer_clustering.model import InitialConfig, initialise
from polymer_clustering.movement import MovementConfig, MovementEngine, retether, simulate_movement
from polymer_clustering.population import PopulationConfig
from polymer_clustering.rates import RateTree
from polymer_clustering.response import ResponseConfig
from polymer_clustering.output import software_metadata
from test_binding import rates, small_world
import test_progress


class OptimisationChecks(unittest.TestCase):
    def test_completed_v010_sources_remain_reusable_only_through_exact_bridges(self):
        current = software_metadata()
        self.assertEqual(current["source_sha256"]["movement.py"], OPTIMISED_MOVEMENT_SHA256)
        self.assertEqual(current["source_sha256"]["rates.py"], OPTIMISED_RATES_SHA256)
        old = copy.deepcopy(current)
        old["source_sha256"].update({"movement.py": PROGRESS_MOVEMENT_SHA256,
                                     "rates.py": PRE_OPTIMISATION_RATES_SHA256})
        for cutoff in (10, None, 20):
            self.assertTrue(numerical_sources_match(old, current,
                                                    {"ps_crowding_threshold": cutoff}))
        for name in ("movement.py", "rates.py", "response.py"):
            edited = copy.deepcopy(old)
            edited["source_sha256"][name] = "unknown predecessor"
            self.assertFalse(numerical_sources_match(edited, current, {}))
            edited = copy.deepcopy(current)
            edited["source_sha256"][name] = "unknown target"
            self.assertFalse(numerical_sources_match(old, edited, {}))
        old["numpy"] = "other version"
        self.assertFalse(numerical_sources_match(old, current, {}))

    def assert_cache(self, engine, exact_trees=True):
        """The readable full rebuild is an independent dependency-update oracle."""
        reference = MovementEngine(copy.deepcopy(engine.state), engine.binding, engine.movement)
        self.assertEqual(engine.groups, reference.groups)
        self.assertEqual(engine.pb_watchers, reference.pb_watchers)
        self.assertEqual(engine.pb_targets, reference.pb_targets)
        self.assertEqual(engine.polymer_bins, reference.polymer_bins)
        self.assertEqual(engine.bs_channels, reference.bs_channels)
        self.assertEqual(engine.polymer_channels, reference.polymer_channels)
        np.testing.assert_array_equal(engine.signal_rates, reference.signal_rates)
        np.testing.assert_array_equal(engine.polymer_rates, reference.polymer_rates)
        np.testing.assert_array_equal(engine.cluster_rates, reference.cluster_rates)
        for actual, expected in zip(engine.trees, reference.trees):
            np.testing.assert_array_equal(actual.values, expected.values)
            self.assertEqual(actual.positive, expected.positive)
            if exact_trees:
                np.testing.assert_array_equal(actual.tree, expected.tree)

    def pb(self, engine, polymer, kind, contact=None):
        """Choose an interior draw for a specified positive PB channel."""
        channels = engine.polymer_channels[polymer]
        index = next(i for i, c in enumerate(channels)
                     if c.kind == kind and (contact is None or c.contact == contact))
        cumulative = np.cumsum([c.propensity for c in channels])
        before = cumulative[index - 1] if index else 0.
        class Draw:
            def random(self):
                return float((before + cumulative[index]) / 2 / cumulative[-1])
        record = engine.event(1, polymer, Draw())
        self.assertEqual(record["reaction"], kind)
        engine.state.validate(20, len(engine.state.polymer_sites))
        self.assert_cache(engine)

    def test_ordered_bulk_sums_are_bit_identical_to_old_sequential_updates(self):
        rng = np.random.default_rng(417)
        for size in (0, 1, 7, 128, 2500, 7000):
            for shuffled in (False, True):
                with self.subTest(size=size, shuffled=shuffled):
                    values = 10.**rng.uniform(-12, 3, size)
                    values[rng.random(size) < .3] = 0
                    order = rng.permutation(size) if shuffled else np.arange(size)
                    reference = RateTree(np.zeros(size))
                    for index in order:
                        reference.update(int(index), float(values[index]))
                    actual = RateTree.from_sequential_updates(values, order if shuffled else None)
                    np.testing.assert_array_equal(reference.values, actual.values)
                    np.testing.assert_array_equal(reference.tree, actual.tree)
                    self.assertEqual(reference.total, actual.total)
                    self.assertEqual(reference.positive, actual.positive)
                    if actual.total:
                        for fraction in (0., .1, .5, .999999):
                            self.assertEqual(reference.select(fraction*reference.total),
                                             actual.select(fraction*actual.total))
        for order in ([0, 0], [0, 2], [.0, 1.]):
            with self.assertRaises(ValueError):
                RateTree.from_sequential_updates(np.ones(2), order)
        for values in ([1, -1], [float("inf")], [[1, 2]]):
            with self.assertRaises(ValueError):
                RateTree.from_sequential_updates(np.array(values))

    def test_unchanged_rate_skips_ancestor_updates_without_changing_bits(self):
        tree = RateTree(np.array([0., .1, .2, 0.]))
        original = tree.tree.copy()
        for _ in range(100):
            tree.update(1, .1)
            tree.update(3, -0.)
        np.testing.assert_array_equal(tree.tree, original)
        self.assertEqual(tree.positive, 2)
        self.assertTrue(np.signbit(tree.values[3]))

    def test_full_site_release_restores_free_and_neighbour_host_channels(self):
        state = BindingState.from_initial(small_world(polymers=33))
        state.polymer_sites[:] = 1
        state.polymer_contacts[:30, 0] = 1
        state.polymer_sites[31] = 2
        state.polymer_contacts[31, 1] = 2
        state.bacteria_polymers.ravel()[1] = 30
        state.bacteria_polymers.ravel()[2] = 1
        state.validate(20, 33)
        engine = MovementEngine(state, rates(), MovementConfig())
        self.assertIn(31, engine.pb_watchers[1])
        self.assertFalse(any(c.kind == "pb_bind" and c.contact == 1
                             for c in engine.polymer_channels[31]))
        self.pb(engine, 0, "pb_unbind", 0)
        for polymer in (30, 31, 32):
            self.assertTrue(any(c.kind == "pb_bind" and c.contact == 1
                                for c in engine.polymer_channels[polymer]))
        self.pb(engine, 30, "pb_bind", 1)
        self.assertFalse(any(c.kind == "pb_bind" and c.contact == 1
                             for c in engine.polymer_channels[31]))

    def test_bridge_release_reanchors_and_crosses_both_ps_cutoff_boundaries(self):
        for cutoff in (10, None, 20):
            with self.subTest(cutoff=cutoff):
                state = BindingState.from_initial(small_world(polymers=19))
                state.polymer_sites[10:] = 2
                state.free_signals.ravel()[0] = 10
                state.free_signals.ravel()[2] = 10
                engine = MovementEngine(state, rates(ps_crowding_threshold=cutoff), MovementConfig())
                self.pb(engine, 0, "pb_bind", 0)
                self.pb(engine, 0, "pb_bind", 2)
                self.assertEqual(engine.groups, ((0, 1),))
                self.pb(engine, 0, "pb_unbind", 0)
                self.assertEqual(engine.groups, ((0,), (1,)))
                self.assertEqual(state.polymer_sites[0], 2)
                self.assertEqual(len(engine.polymer_bins[0]), 9)
                self.assertEqual(len(engine.polymer_bins[2]), 10)
                if cutoff == 10:
                    self.assertTrue(any(c.kind == "ps_bind" for c in engine.polymer_channels[1]))
                    self.assertFalse(any(c.kind == "ps_bind" for c in engine.polymer_channels[10]))
                self.pb(engine, 0, "pb_unbind", 1)
                self.assertEqual(state.polymer_sites[0], 2)

    def test_redundant_bridge_release_retains_the_connected_component(self):
        state = BindingState.from_initial(small_world(polymers=2))
        engine = MovementEngine(state, rates(), MovementConfig())
        for polymer in (0, 1):
            self.pb(engine, polymer, "pb_bind", 0)
            self.pb(engine, polymer, "pb_bind", 2)
        self.pb(engine, 0, "pb_unbind", 0)
        self.assertEqual(engine.groups, ((0, 1),))
        self.pb(engine, 1, "pb_unbind", 1)
        self.assertEqual(engine.groups, ((0,), (1,)))

    def test_polymer_hop_updates_dependencies_and_lazy_insertion_order(self):
        state = BindingState.from_initial(small_world(polymers=3))
        engine = MovementEngine(state, rates(), MovementConfig())
        class East:
            def random(self):
                return .75  # Equal positive +y/+x rates at the upper-left site.
        engine.event(3, 0, East())
        self.assertEqual(state.polymer_sites[0], 1)
        self.assertEqual(engine.pb_targets[0], (1,))
        self.assertNotIn(0, engine.pb_watchers[0])
        self.assertIn(0, engine.pb_watchers[1])
        self.assert_cache(engine, exact_trees=False)
        self.pb(engine, 1, "pb_bind", 0)
        self.assertFalse(engine.chemical_order_dirty)
        np.testing.assert_array_equal(engine.chemical_order, [1, 2, 0])

    def test_mixed_growth_production_and_transport_replay_full_rebuild_oracle(self):
        def legacy_pb(engine, polymer, site, kind):
            retether(engine.state)
            engine.rebuild()
        checker = test_progress.ProgressChecks()
        seen = set()
        for seed in (21, 22):
            initial = initialise(InitialConfig(grid_size=8, bacteria=4, polymers=160,
                                               strain="BB170", seed=seed,
                                               rng_scheme="independent_v1"))
            config = rates(duration=10, sample_interval=2, k_pb_on=.2, k_pb_off=.1,
                           k_ps_on=.1, k_ps_off=.1)
            movement = MovementConfig(d_s=.2, d_p=.2, d_b=.2, spatial_frames=3,
                                      event_log="all", validate_every=17)
            population = PopulationConfig(doubling_time=3, division_width=0,
                                          population_cap=12, production_rate=4)
            def run():
                return simulate_movement(initial, config, movement, seed,
                    response=ResponseConfig(), strain="BB170", population=population)
            with patch.object(MovementEngine, "refresh_after_pb", legacy_pb):
                reference = run()
            actual = run()
            checker.assert_same(reference, actual)
            seen.update(actual.event_counts)
        self.assertTrue({"division", "signal_production", "pb_bind", "pb_unbind",
                         "polymer_hop", "cluster_move"} <= seen, seen)


if __name__ == "__main__":
    unittest.main()
