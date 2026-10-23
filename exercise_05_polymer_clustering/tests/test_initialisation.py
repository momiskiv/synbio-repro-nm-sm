# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Scientific invariants and reproducibility checks for initialisation."""

import unittest

import numpy as np

from polymer_clustering.model import InitialConfig, initialise


class InitialisationTests(unittest.TestCase):
    def test_default_footprints_and_counts(self) -> None:
        # Every cell must occupy four sites bearing its own ID, without overlap.
        state = initialise(InitialConfig())
        self.assertEqual(state.occupancy.shape, (50, 50))
        self.assertEqual(len(state.bacteria), 25)
        self.assertEqual(np.count_nonzero(state.occupancy >= 0), 100)
        for bacterium in state.bacteria:
            self.assertTrue(0 <= bacterium.x <= 48)
            self.assertTrue(0 <= bacterium.y <= 48)
            footprint = state.occupancy[bacterium.y : bacterium.y + 2,
                                        bacterium.x : bacterium.x + 2]
            np.testing.assert_array_equal(footprint, np.full((2, 2), bacterium.identifier))
            self.assertEqual(np.count_nonzero(state.occupancy == bacterium.identifier), 4)
        self.assertEqual(int(state.free_signals.sum()), 10_000)
        self.assertEqual(int(state.free_polymers.sum()), 0)
        self.assertTrue(np.all(state.free_signals >= 0))

    def test_seed_repeats_state(self) -> None:
        # Reusing a seed must reproduce the full state, not just its totals.
        first = initialise(InitialConfig(polymers=1000))
        second = initialise(InitialConfig(polymers=1000))
        self.assertEqual(first.bacteria, second.bacteria)
        for name in ("occupancy", "free_signals", "free_polymers"):
            np.testing.assert_array_equal(getattr(first, name), getattr(second, name))
        different = initialise(InitialConfig(seed=2027))
        self.assertNotEqual(first.bacteria, different.bacteria)

    def test_bb170_and_polymer_counts(self) -> None:
        # Check strain-specific starting signals and exact polymer allocation.
        state = initialise(InitialConfig(strain="BB170", polymers=7000))
        self.assertEqual(int(state.free_signals.sum()), 0)
        self.assertEqual(int(state.free_polymers.sum()), 7000)

    def test_empty_and_smallest_states(self) -> None:
        # Boundary cases exercise both the empty grid and a single valid anchor.
        empty = initialise(InitialConfig(grid_size=2, bacteria=0))
        self.assertTrue(np.all(empty.occupancy == -1))
        full = initialise(InitialConfig(grid_size=2, bacteria=1))
        np.testing.assert_array_equal(full.occupancy, np.zeros((2, 2)))

    def test_invalid_inputs(self) -> None:
        # Reject impossible geometry, unsupported strains, and invalid types.
        for settings in (
            {"grid_size": 1}, {"bacteria": -1}, {"bacteria": 81},
            {"grid_size": 3, "bacteria": 2}, {"polymers": -1},
            {"strain": "unknown"}, {"seed": -1}, {"seed": 2**64},
            {"bacteria": True}, {"polymers": 1.5},
        ):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                InitialConfig(**settings)


if __name__ == "__main__":
    unittest.main()

