# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Analytical response solutions, sensing clocks, and passive coupling checks."""

import unittest

import numpy as np

from polymer_clustering.binding import BindingState, apply_reaction
from polymer_clustering.movement import MovementConfig, move_cluster, simulate_movement
from polymer_clustering.response import (
    ResponseConfig, ResponseTracker, advance_expression, light_from_expression, response_summary,
)
from test_binding import channel, rates, small_world
from test_movement import empty_world


class ResponseTests(unittest.TestCase):
    def test_synthesis_matches_closed_form_solution(self) -> None:
        config = ResponseConfig()
        value = advance_expression(np.array([0.0]), np.array([38]), np.array([20.0]), config)
        self.assertAlmostEqual(value[0], 38 * (1 - np.exp(-1)))
        # Two intervals with constant occupancy must equal one combined interval.
        split = advance_expression(value, np.array([38]), np.array([10.0]), config)
        combined = advance_expression(np.array([0.0]), np.array([38]), np.array([30.0]), config)
        np.testing.assert_allclose(split, combined, rtol=1e-14)

    def test_decay_sign_crossing_and_zero_sensing(self) -> None:
        config = ResponseConfig()
        early = advance_expression(np.array([50.0]), np.array([10]), np.array([10.0]), config)
        self.assertAlmostEqual(early[0], 50 * np.exp(-0.1))
        late = advance_expression(np.array([50.0]), np.array([10]), np.array([1000.0]), config)
        self.assertEqual(late[0], 10)  # The declining branch stops at the sensed level.
        zero = advance_expression(np.array([50.0]), np.array([0]), np.array([100.0]), config)
        self.assertAlmostEqual(zero[0], 50 * np.exp(-1))
        still = advance_expression(np.array([20.0]), np.array([20]), np.array([100.0]), config)
        self.assertEqual(still[0], 20)

    def test_zero_rates_and_tiny_intervals(self) -> None:
        config = ResponseConfig(k_syn=0, k_deg=0)
        values = np.array([0.0, 50.0])
        np.testing.assert_array_equal(
            advance_expression(values, np.array([20, 0]), np.array([100.0, 100.0]), config),
            values,
        )
        tiny = advance_expression(np.array([0.0]), np.array([1]), np.array([1e-12]), ResponseConfig())
        self.assertGreater(tiny[0], 0)
        self.assertAlmostEqual(tiny[0] / 1e-12, 0.05, places=12)

    def test_light_midpoint_baseline_and_asymptotes(self) -> None:
        config = ResponseConfig()
        values = light_from_expression(np.array([-1e6, 38.0, 1e6]), 38, 2500)
        np.testing.assert_allclose(values, [0, 2500, 5000])
        self.assertTrue(np.all(np.isfinite(values)))
        zero_expression = light_from_expression(np.array([0.0]), 38, 2500)[0]
        self.assertAlmostEqual(zero_expression, 5000 / (1 + np.exp(9.5)))
        self.assertGreater(zero_expression, 0)
        self.assertEqual(config.midpoint("MM32"), 38)
        self.assertEqual(config.midpoint("BB170"), 15)
        self.assertEqual(ResponseConfig(response_midpoint=22).midpoint("MM32"), 22)

    def test_tracker_updates_old_occupancy_before_binding_change_and_moves_with_cell(self) -> None:
        state = BindingState.from_initial(small_world(polymers=0))
        tracker = ResponseTracker(state, ResponseConfig(), "MM32")
        apply_reaction(state, channel(state, "bs_bind"))
        tracker.change_sensing(0, 1, 10)
        self.assertEqual(tracker.observe(10)["expression"][0], 0)
        move_cluster(state, (0,), 1)
        tracker.validate(state, 20)
        before_release = 1 - np.exp(-0.5)
        self.assertAlmostEqual(tracker.observe(20)["expression"][0], before_release)
        apply_reaction(state, channel(state, "bs_unbind"))
        tracker.change_sensing(0, -1, 20)
        self.assertAlmostEqual(tracker.observe(30)["expression"][0], before_release * np.exp(-0.1))
        tracker.validate(state, 30)

    def test_population_light_is_sum_of_cell_responses(self) -> None:
        expression = np.array([0.0, 20.0])
        lights = light_from_expression(expression, 38, 2500)
        values = {"bacterium": np.array([0, 1]), "sensed": np.array([0, 20]),
                  "expression": expression, "luminescence": lights}
        summary = response_summary(values, 0)
        self.assertAlmostEqual(summary["luminescence_total"], float(lights.sum()))
        # Nonlinearity means f(mean expression) generally differs from mean f(E).
        at_mean = light_from_expression(np.array([10.0]), 38, 2500)[0]
        self.assertNotAlmostEqual(summary["luminescence_mean"], at_mean)

    def test_response_is_passive_and_does_not_change_stochastic_history(self) -> None:
        initial = small_world(polymers=2)
        binding = rates(duration=10, sample_interval=1)
        movement = MovementConfig(spatial_frames=3, event_log="all")
        plain = simulate_movement(initial, binding, movement, 4)
        responding = simulate_movement(initial, binding, movement, 4, response=ResponseConfig())
        self.assertEqual(plain.events, responding.events)
        self.assertEqual(plain.trajectory, responding.trajectory)
        self.assertEqual(plain.state.bacteria, responding.state.bacteria)
        np.testing.assert_array_equal(plain.state.free_signals, responding.state.free_signals)
        for row, counts in zip(responding.response_trajectory, responding.trajectory):
            self.assertEqual(row["sensed_total"], counts["bacteria_bound_signals"])

    def test_sampling_frames_and_logs_do_not_change_expression(self) -> None:
        initial = small_world(polymers=2)
        first = simulate_movement(initial, rates(duration=10, sample_interval=1),
                                  MovementConfig(spatial_frames=3, event_log="all"),
                                  4, response=ResponseConfig())
        second = simulate_movement(initial, rates(duration=10, sample_interval=0.3),
                                   MovementConfig(spatial_frames=8, event_log="compact"),
                                   4, response=ResponseConfig())
        np.testing.assert_array_equal(first.response_final["expression"],
                                      second.response_final["expression"])
        np.testing.assert_array_equal(first.response_final["luminescence"],
                                      second.response_final["luminescence"])
        self.assertEqual(first.response_trajectory[-1], second.response_trajectory[-1])

    def test_empty_population_and_no_signal_baseline(self) -> None:
        empty = simulate_movement(empty_world(signals=0), rates(duration=10),
                                   MovementConfig(spatial_frames=2), 1, response=ResponseConfig())
        self.assertEqual(empty.response_trajectory[-1]["luminescence_total"], 0)
        self.assertEqual(empty.response_trajectory[-1]["luminescence_mean"], 0)
        self.assertEqual(empty.response_final["expression"].size, 0)
        no_signal = simulate_movement(small_world(polymers=0, signals=0), rates(duration=10),
                                       MovementConfig(spatial_frames=2), 1, response=ResponseConfig())
        np.testing.assert_array_equal(no_signal.response_final["expression"], [0, 0])
        self.assertGreater(no_signal.response_trajectory[-1]["luminescence_total"], 0)

    def test_event_cap_response_stops_at_actual_time(self) -> None:
        result = simulate_movement(small_world(), rates(duration=10, max_events=1),
                                   MovementConfig(spatial_frames=2), 2, response=ResponseConfig())
        self.assertEqual(result.status, "event_limit")
        self.assertLess(result.time, 10)
        self.assertEqual(result.response_trajectory[-1]["time"], result.time)
        self.assertTrue(all(row["time"] <= result.time for row in result.cell_response))
        np.testing.assert_array_equal(result.frames[-1].expression,
                                      result.response_final["expression"])

    def test_invalid_response_values(self) -> None:
        for settings in ({"k_deg": -0.01}, {"k_syn": float("nan")}, {"light_scale": 0},
                         {"response_initial": 61}, {"response_midpoint": True},
                         {"light_map_max": 0}):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                ResponseConfig(**settings)
        with self.assertRaises(ValueError):
            advance_expression(np.array([0.0]), np.array([1]), np.array([-1.0]), ResponseConfig())


if __name__ == "__main__":
    unittest.main()
