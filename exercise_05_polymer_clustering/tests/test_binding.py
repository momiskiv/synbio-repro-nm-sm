# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Reaction kinetics, conservation, contact geometry, and stochastic timing checks."""

import unittest

import numpy as np

from polymer_clustering.binding import (
    BindingConfig, BindingState, apply_reaction, neighbour_contacts,
    reaction_channels, sample_counts, simulate_binding,
)
from polymer_clustering.model import Bacterium, InitialConfig, InitialState, initialise


def small_world(polymers: int = 1, signals: int = 20, diagonal: bool = False) -> InitialState:
    """Build a controlled two-cell world to test contacts and exact reaction rates."""
    occupancy = np.full((4, 4), -1, dtype=np.int32)
    second_y = 2 if diagonal else 0
    bacteria = (Bacterium(0, 0, 0), Bacterium(1, 2, second_y))
    for bacterium in bacteria:
        occupancy[bacterium.y:bacterium.y + 2, bacterium.x:bacterium.x + 2] = bacterium.identifier
    free_signals = np.zeros((4, 4), dtype=np.int64)
    free_signals[0, 0] = signals
    free_polymers = np.zeros_like(free_signals)
    free_polymers[0, 0] = polymers
    return InitialState(occupancy, bacteria, free_signals, free_polymers)


def rates(**overrides: float) -> BindingConfig:
    """Explicit test kinetics; the polymer values are not inferred from the paper."""
    return BindingConfig(**{
        "k_pb_on": 0.005, "k_pb_off": 0.001,
        "k_ps_on": 0.01, "k_ps_off": 0.001, **overrides,
    })


def channel(state: BindingState, kind: str, config: BindingConfig | None = None):
    """Select one currently possible event for an explicit reaction sequence."""
    return next(item for item in reaction_channels(
        state, config or rates(), neighbour_contacts(state)
    ) if item.kind == kind)


class BindingTests(unittest.TestCase):
    def test_published_bs_mass_action_and_reversal(self) -> None:
        state = BindingState.from_initial(small_world(polymers=0, signals=4))
        forward = channel(state, "bs_bind")
        self.assertAlmostEqual(forward.propensity, 0.001 * 15 * 4)
        apply_reaction(state, forward)
        self.assertEqual(state.free_signals[0, 0], 3)
        self.assertEqual(state.bacteria_signals[0, 0], 1)
        backward = channel(state, "bs_unbind")
        self.assertAlmostEqual(backward.propensity, 0.002)
        apply_reaction(state, backward)
        state.validate(4, 0)
        self.assertEqual(state.free_signals[0, 0], 4)

    def test_ps_saturation_reversal_and_conservation(self) -> None:
        state = BindingState.from_initial(small_world())
        self.assertAlmostEqual(channel(state, "ps_bind").propensity, 0.01 * 5 * 20)
        for _ in range(5):
            apply_reaction(state, channel(state, "ps_bind"))
            state.validate(20, 1)
        self.assertFalse(any(r.kind == "ps_bind" for r in reaction_channels(
            state, rates(), neighbour_contacts(state))))
        self.assertAlmostEqual(channel(state, "ps_unbind").propensity, 0.001 * 5)
        for _ in range(5):
            apply_reaction(state, channel(state, "ps_unbind"))
        state.validate(20, 1)
        self.assertEqual(state.free_signals[0, 0], 20)

    def test_bs_site_capacity_and_locality(self) -> None:
        state = BindingState.from_initial(small_world(polymers=0))
        for _ in range(15):
            apply_reaction(state, channel(state, "bs_bind"))
        self.assertEqual(state.bacteria_signals.sum(), 15)
        self.assertFalse(any(r.kind == "bs_bind" for r in reaction_channels(
            state, rates(), neighbour_contacts(state))))
        # Signals outside bacterial footprints cannot bind a bacterial receptor.
        state.free_signals[0, 0] -= 1
        state.free_signals[3, 0] += 1
        state.validate(20, 0)
        self.assertFalse(any(r.kind == "bs_bind" for r in reaction_channels(
            state, rates(), neighbour_contacts(state))))

    def test_bridge_identity_and_both_end_release(self) -> None:
        state = BindingState.from_initial(small_world())
        first = channel(state, "pb_bind")
        self.assertAlmostEqual(first.propensity, 0.005 * 2 * 30)
        apply_reaction(state, first)
        possibilities = [r for r in reaction_channels(state, rates(), neighbour_contacts(state))
                         if r.kind == "pb_bind"]
        self.assertEqual({r.contact for r in possibilities}, {2, 6})
        self.assertTrue(all(r.propensity == 0.005 * 30 for r in possibilities))
        apply_reaction(state, possibilities[0])
        state.validate(20, 1)
        counts = sample_counts(state, 0)
        self.assertEqual(counts["bridges"], 1)
        self.assertEqual(counts["largest_cluster"], 2)
        self.assertEqual(counts["pb_bonds"], 2)
        self.assertEqual(counts["attached_polymers"], 1)
        self.assertFalse(any(r.kind == "pb_bind" for r in reaction_channels(
            state, rates(), neighbour_contacts(state))))
        # Releasing either end breaks the bridge; releasing the last returns the
        # polymer to the unattached population without destroying its identity.
        for _ in range(2):
            apply_reaction(state, channel(state, "pb_unbind"))
            state.validate(20, 1)
        self.assertEqual(sample_counts(state, 0)["free_polymers"], 1)

    def test_diagonal_bacteria_cannot_bridge(self) -> None:
        state = BindingState.from_initial(small_world(diagonal=True))
        apply_reaction(state, channel(state, "pb_bind"))
        self.assertFalse(any(r.kind == "pb_bind" for r in reaction_channels(
            state, rates(), neighbour_contacts(state))))

    def test_crowding_cutoff_preserves_bound_signals_and_release(self) -> None:
        nine = BindingState.from_initial(small_world(polymers=9))
        self.assertTrue(any(r.kind == "ps_bind" for r in reaction_channels(
            nine, rates(), neighbour_contacts(nine))))
        ten = BindingState.from_initial(small_world(polymers=10))
        ten.free_signals[0, 0] -= 1
        ten.polymer_signals[0] = 1
        self.assertFalse(any(r.kind == "ps_bind" for r in reaction_channels(
            ten, rates(), neighbour_contacts(ten))))
        apply_reaction(ten, channel(ten, "ps_unbind"))
        ten.validate(20, 10)
        self.assertEqual(ten.free_signals[0, 0], 20)

    def test_pb_site_capacity_and_ps_independence(self) -> None:
        state = BindingState.from_initial(small_world(polymers=40))
        # Fill the 30 polymer sites at one bacterial lattice site. Force only
        # first attachments here; neighbouring contact sites remain distinct.
        for polymer in range(30):
            chosen = next(r for r in reaction_channels(state, rates(), neighbour_contacts(state))
                          if r.kind == "pb_bind" and r.target == polymer and r.contact == 0)
            apply_reaction(state, chosen)
        state.validate(20, 40)
        self.assertFalse(any(r.kind == "pb_bind" and r.contact == 0 for r in reaction_channels(
            state, rates(), neighbour_contacts(state))))
        # Binding a signal to one polymer does not consume its bacterial sites.
        other = BindingState.from_initial(small_world())
        apply_reaction(other, channel(other, "ps_bind"))
        apply_reaction(other, channel(other, "pb_bind"))
        apply_reaction(other, channel(other, "pb_unbind"))
        other.validate(20, 1)
        self.assertEqual(other.polymer_signals[0], 1)

    def test_seed_replay_conservation_and_initial_state_preservation(self) -> None:
        initial = small_world(polymers=3, signals=50)
        config = rates(duration=20, sample_interval=0.7)
        first = simulate_binding(initial, config, 5)
        second = simulate_binding(initial, config, 5)
        self.assertEqual(first.status, "complete")
        self.assertEqual(first.events, second.events)
        self.assertEqual(first.trajectory, second.trajectory)
        self.assertEqual(first.trajectory[-1]["time"], 20)
        self.assertTrue(all(row["total_signals"] == 50 for row in first.trajectory))
        self.assertTrue(all(row["free_polymers"] + row["attached_polymers"] == 3
                            for row in first.trajectory))
        self.assertEqual(initial.free_signals[0, 0], 50)
        self.assertEqual(initial.free_polymers[0, 0], 3)

    def test_sampling_interval_does_not_change_event_history(self) -> None:
        initial = small_world()
        first = simulate_binding(initial, rates(duration=10, sample_interval=1), 7)
        second = simulate_binding(initial, rates(duration=10, sample_interval=0.3), 7)
        self.assertEqual(first.events, second.events)

    def test_no_event_after_horizon_and_explicit_event_cap(self) -> None:
        initial = small_world()
        short = simulate_binding(initial, rates(duration=1e-15), 1)
        self.assertEqual(short.events, [])
        self.assertEqual(short.time, 1e-15)
        capped = simulate_binding(initial, rates(duration=100, max_events=1), 1)
        self.assertEqual(capped.status, "event_limit")
        self.assertEqual(len(capped.events), 1)
        self.assertEqual(capped.trajectory[-1]["time"], capped.time)
        self.assertLess(capped.time, 100)

    def test_zero_propensity_bb170_and_disabled_channels(self) -> None:
        initial = initialise(InitialConfig(strain="BB170"))
        result = simulate_binding(initial, BindingConfig(duration=10), 1)
        self.assertEqual(result.events, [])
        self.assertEqual(result.time, 10)
        self.assertTrue(all(row["total_signals"] == 0 for row in result.trajectory))
        result = simulate_binding(small_world(), rates(
            k_bs_on=0, k_bs_off=0, k_pb_on=0, k_pb_off=0, k_ps_on=0, k_ps_off=0
        ), 1)
        self.assertEqual(result.events, [])

    def test_invalid_rates_and_missing_polymer_kinetics(self) -> None:
        for settings in ({"k_bs_on": -1}, {"duration": 0}, {"sample_interval": 0},
                         {"k_ps_on": float("nan")}, {"k_pb_off": float("inf")},
                         {"max_events": True}, {"k_bs_off": True},
                         {"duration": 100001, "sample_interval": 1}):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                BindingConfig(**settings)
        with self.assertRaises(ValueError):
            simulate_binding(small_world(), BindingConfig(), 1)

    def test_invariant_checker_rejects_lost_signals_and_bad_contacts(self) -> None:
        state = BindingState.from_initial(small_world())
        state.free_signals[0, 0] -= 1
        with self.assertRaises(ValueError):
            state.validate(20, 1)
        state.free_signals[0, 0] += 1
        state.polymer_contacts[0, 0] = 15  # No bacterium occupies this site.
        with self.assertRaises(ValueError):
            state.validate(20, 1)

    def test_competing_reactions_and_exponential_waiting_times(self) -> None:
        # One signal competes between 15 bacterial sites and 5 polymer sites.
        # Rates 0.1 and 0.2 give total 2.5, P(BS first)=0.6, mean wait=0.4.
        # Fixed seeds make this statistical check reproducible. The tolerances
        # cover sampling variation while rejecting substantial kinetic errors.
        config = rates(k_bs_on=0.1, k_bs_off=0, k_ps_on=0.2,
                       k_ps_off=0, k_pb_on=0, k_pb_off=0,
                       duration=100, max_events=1)
        kinds, waits = [], []
        for seed in range(600):
            result = simulate_binding(small_world(signals=1), config, seed)
            kinds.append(result.events[0]["reaction"])
            waits.append(result.events[0]["time"])
        self.assertAlmostEqual(kinds.count("bs_bind") / len(kinds), 0.6, delta=0.06)
        self.assertAlmostEqual(float(np.mean(waits)), 0.4, delta=0.06)

    def test_reversible_one_signal_matches_analytic_solution(self) -> None:
        # One signal is either free or bound. The free-to-bound rate is
        # 15 * 0.1 = 1.5 and the reverse rate is 0.5. Starting free gives
        # P(bound at t) = 0.75 * (1 - exp(-2*t)), hence approximately 0.75 at t=10.
        config = BindingConfig(k_bs_on=0.1, k_bs_off=0.5, duration=10,
                               sample_interval=10)
        bound = []
        for seed in range(600):
            result = simulate_binding(small_world(polymers=0, signals=1), config, seed)
            bound.append(int(result.state.bacteria_signals.sum()))
        expected = 0.75 * (1 - np.exp(-20))
        self.assertAlmostEqual(float(np.mean(bound)), expected, delta=0.06)


if __name__ == "__main__":
    unittest.main()
