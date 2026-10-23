# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Division appointments, inheritance, Poisson production, and conservation checks."""

import unittest

import numpy as np

from polymer_clustering.binding import BindingConfig, BindingState, apply_reaction
from polymer_clustering.model import InitialConfig, initialise
from polymer_clustering.movement import MovementConfig, MovementEngine, simulate_movement
from polymer_clustering.population import (
    DivisionSchedule, PopulationConfig, nearest_empty_anchor, produce_signal,
)
from polymer_clustering.response import ResponseConfig, ResponseTracker
from test_binding import channel, rates, small_world


def quiet_movement(**overrides):
    """Isolate population rules from transport without changing their event clock."""
    return MovementConfig(**{"d_s": 0, "d_p": 0, "d_b": 0, "spatial_frames": 3, **overrides})


def no_binding(duration=6, **overrides):
    return BindingConfig(**{"duration": duration, "sample_interval": 1,
                            "k_bs_on": 0, "k_bs_off": 0, **overrides})


class PopulationTests(unittest.TestCase):
    def test_triangular_lifetimes_match_source_mean_and_full_width_choice(self):
        config = PopulationConfig()
        rng = np.random.default_rng(32)
        lifetimes = np.array([config.lifetime(rng) for _ in range(20000)])
        self.assertGreaterEqual(lifetimes.min(), 5000)
        self.assertLessEqual(lifetimes.max(), 25000)
        self.assertAlmostEqual(lifetimes.mean(), 15000, delta=100)
        # A symmetric triangle's variance is width squared / 24, not the
        # variance of a uniform distribution on the same interval.
        self.assertAlmostEqual(lifetimes.var(), 20000**2 / 24, delta=600000)
        self.assertEqual(PopulationConfig(doubling_time=2, division_width=0).lifetime(rng), 2)

    def test_nearest_placement_retains_parent_and_all_attachment_pools(self):
        state = BindingState.from_initial(small_world())
        for kind in ("bs_bind", "pb_bind", "ps_bind"):
            apply_reaction(state, channel(state, kind))
        before = [a.copy() for a in (state.bacteria_signals, state.bacteria_polymers,
                                     state.polymer_signals, state.polymer_contacts,
                                     state.free_signals)]
        config = PopulationConfig(doubling_time=2, division_width=0, population_cap=3)
        schedule = DivisionSchedule(state, config, np.random.default_rng(1))
        record = schedule.divide(state, 2, np.random.default_rng(2))
        daughter = state.bacteria[-1]
        self.assertEqual(record["reaction"], "division")
        self.assertEqual((daughter.x, daughter.y), (0, 2))  # unique nearest anchor to cell 0
        self.assertEqual(schedule.next_time, np.inf)
        for old, actual in zip(before, (state.bacteria_signals, state.bacteria_polymers,
                                        state.polymer_signals, state.polymer_contacts,
                                        state.free_signals)):
            np.testing.assert_array_equal(old, actual)
        state.validate(20, 1)
        # A complete rebuild must include the daughter's newly available receptors.
        engine = MovementEngine(state, rates(), quiet_movement())
        self.assertEqual(len(engine.groups), 3)
        self.assertEqual(schedule.lineage[-1]["parent"], 0)

    def test_nearest_ties_are_reachable_and_blocked_placement_retries(self):
        state = BindingState.from_initial(initialise(InitialConfig(grid_size=6, bacteria=1,
                                                                   strain="BB170")))
        # Move the only cell to an interior anchor: all four distance-2 choices exist.
        from polymer_clustering.model import Bacterium
        state.occupancy[:] = -1
        state.bacteria = (Bacterium(0, 2, 2),)
        state.occupancy[2:4, 2:4] = 0
        rng = np.random.default_rng(4)
        anchors = {nearest_empty_anchor(state, state.bacteria[0], rng) for _ in range(100)}
        self.assertEqual(anchors, {(0, 2), (4, 2), (2, 0), (2, 4)})
        jammed = BindingState.from_initial(initialise(InitialConfig(grid_size=2, bacteria=1,
                                                                    strain="BB170")))
        schedule = DivisionSchedule(jammed, PopulationConfig(doubling_time=2, division_width=0,
                                                             division_retry=3), rng)
        self.assertEqual(schedule.divide(jammed, 2, rng)["reaction"], "division_blocked")
        self.assertEqual(schedule.next_time, 5)
        self.assertEqual(len(jammed.bacteria), 1)

    def test_age_scheduling_simultaneous_divisions_and_cap(self):
        initial = initialise(InitialConfig(grid_size=8, bacteria=1, strain="BB170"))
        result = simulate_movement(initial, no_binding(duration=6), quiet_movement(), 7,
                                   ResponseConfig(), "BB170",
                                   PopulationConfig(doubling_time=2, division_width=0,
                                                    population_cap=4))
        self.assertEqual(result.event_counts["division"], 3)
        self.assertEqual([row["time"] for row in result.events], [2, 4, 4])
        self.assertEqual([row["bacteria"] for row in result.population_trajectory],
                         [1, 1, 2, 2, 4, 4, 4])
        self.assertEqual(result.time, 6)
        self.assertEqual(len(result.response_final["bacterium"]), 4)
        result.state.validate(result.produced_signals, 0)

    def test_daughter_response_clock_starts_at_birth_and_parent_is_retained(self):
        state = BindingState.from_initial(small_world())
        apply_reaction(state, channel(state, "bs_bind"))
        tracker = ResponseTracker(state, ResponseConfig(response_initial=10), "MM32")
        before = tracker.observe(20)
        tracker.add_daughter(2, 20, 8)
        after = tracker.observe(20)
        np.testing.assert_array_equal(before["expression"], after["expression"][:2])
        self.assertEqual(after["sensed"][-1], 0)
        self.assertEqual(after["expression"][-1], 8)
        self.assertAlmostEqual(tracker.observe(21)["expression"][-1], 8 * np.exp(-0.01))

    def test_production_site_and_analytical_poisson_mean(self):
        initial = initialise(InitialConfig(grid_size=6, bacteria=2, strain="BB170"))
        state = BindingState.from_initial(initial)
        rng = np.random.default_rng(4)
        for _ in range(1000):
            record = produce_signal(state, rng)
            self.assertEqual(state.occupancy[record["y"], record["x"]], record["bacterium"])
        self.assertEqual(state.free_signals.sum(), 1000)
        # n=2, rate=0.5 per cell, t=4 gives a Poisson mean AND variance of 4.
        counts = [simulate_movement(initial, no_binding(duration=4), quiet_movement(), seed,
                                   strain="BB170", population=PopulationConfig(
                                       growth=False, production_rate=0.5)).produced_signals
                  for seed in range(300)]
        self.assertAlmostEqual(np.mean(counts), 4, delta=0.4)
        self.assertAlmostEqual(np.var(counts), 4, delta=0.8)
        self.assertEqual(PopulationConfig().signal_rate("MM32"), 0)
        self.assertEqual(PopulationConfig().signal_rate("BB170"), 0.005)

    def test_growth_changes_production_rate_and_every_pool_remains_conserved(self):
        initial = initialise(InitialConfig(grid_size=8, bacteria=1, strain="BB170", polymers=3))
        config = PopulationConfig(doubling_time=2, division_width=0, population_cap=4,
                                  production_rate=2)
        result = simulate_movement(initial, rates(duration=6, sample_interval=0.5),
                                   MovementConfig(d_s=0.1, d_b=0.05, validate_every=1,
                                                  spatial_frames=5), 7, ResponseConfig(),
                                   "BB170", config)
        self.assertGreater(result.produced_signals, 0)
        self.assertEqual(result.event_counts["signal_production"], result.produced_signals)
        self.assertEqual(len(result.state.bacteria), 4)
        self.assertTrue(all(row["expected_total_signals"] == row["actual_total_signals"]
                            for row in result.population_trajectory))
        result.state.validate(result.produced_signals, 3)

    def test_production_hazard_follows_population_exposure(self):
        initial = initialise(InitialConfig(grid_size=8, bacteria=1, strain="BB170"))
        config = PopulationConfig(doubling_time=2, division_width=0, population_cap=4,
                                  production_rate=0.5)
        # N(t)=1 on [0,2), 2 on [2,4), and 4 on [4,6]. Integrated cell-time
        # is 14; Poisson production must therefore have mean 0.5*14=7.
        counts = [simulate_movement(initial, no_binding(duration=6), quiet_movement(), seed,
                                   strain="BB170", population=config).produced_signals
                  for seed in range(300)]
        self.assertAlmostEqual(np.mean(counts), 7, delta=0.6)

    def test_new_receptors_join_chemistry_and_response_after_birth(self):
        initial = small_world(polymers=0, signals=20)
        initial.free_signals[:] = 0
        initial.free_signals[2, 0] = 20  # future daughter's footprint, outside both parents
        result = simulate_movement(initial, no_binding(duration=2, k_bs_on=1),
                                   quiet_movement(), 4, ResponseConfig(), "MM32",
                                   PopulationConfig(doubling_time=1, division_width=0,
                                                    population_cap=3))
        binds = [row for row in result.events if row["reaction"] == "bs_bind"]
        self.assertTrue(binds)
        self.assertTrue(all(row["bacterium"] == 2 and row["time"] > 1 for row in binds))
        self.assertGreater(result.response_final["expression"][-1], 0)
        result.state.validate(20, 0)

    def test_empty_population_has_no_events_or_production(self):
        initial = initialise(InitialConfig(grid_size=2, bacteria=0, strain="BB170"))
        result = simulate_movement(initial, no_binding(), quiet_movement(), 0,
                                   ResponseConfig(), "BB170", PopulationConfig())
        self.assertEqual(result.total_events, 0)
        self.assertEqual(result.produced_signals, 0)
        self.assertEqual(result.lineage, [])
        self.assertEqual(result.population_trajectory[-1]["bacteria"], 0)
        self.assertEqual(result.response_trajectory[-1]["luminescence_total"], 0)

    def test_observations_do_not_change_births_production_or_response(self):
        initial = initialise(InitialConfig(grid_size=8, bacteria=2, strain="BB170", polymers=2))
        config = PopulationConfig(doubling_time=3, division_width=2, population_cap=6,
                                  production_rate=1)
        a = simulate_movement(initial, rates(duration=8, sample_interval=1),
                              quiet_movement(spatial_frames=3), 8, ResponseConfig(), "BB170", config)
        b = simulate_movement(initial, rates(duration=8, sample_interval=0.3),
                              quiet_movement(spatial_frames=11), 8, ResponseConfig(), "BB170", config)
        self.assertEqual(a.events, b.events)
        self.assertEqual(a.lineage, b.lineage)
        for key in a.response_final:
            np.testing.assert_array_equal(a.response_final[key], b.response_final[key])
        np.testing.assert_array_equal(a.state.free_signals, b.state.free_signals)

    def test_disabled_population_preserves_diagnostic_history(self):
        initial = small_world()
        a = simulate_movement(initial, rates(duration=5), MovementConfig(), 9, ResponseConfig())
        b = simulate_movement(initial, rates(duration=5), MovementConfig(), 9, ResponseConfig(),
                              population=PopulationConfig(growth=False, production_rate=0))
        self.assertEqual(a.events, b.events)
        np.testing.assert_array_equal(a.state.free_signals, b.state.free_signals)
        np.testing.assert_array_equal(a.response_final["expression"], b.response_final["expression"])

    def test_event_cap_records_actual_division_time_without_future_frames(self):
        initial = initialise(InitialConfig(grid_size=8, bacteria=1, strain="BB170"))
        result = simulate_movement(initial, no_binding(duration=10, max_events=1),
                                   quiet_movement(), 0, ResponseConfig(), "BB170",
                                   PopulationConfig(doubling_time=2, division_width=0))
        self.assertEqual((result.status, result.time), ("event_limit", 2))
        self.assertEqual(result.frames[-1].time, 2)
        self.assertEqual(result.response_trajectory[-1]["time"], 2)
        self.assertEqual(len(result.state.bacteria), 2)

    def test_invalid_population_controls(self):
        for settings in ({"doubling_time": 0}, {"division_width": 30002},
                         {"division_retry": 0}, {"population_cap": 81},
                         {"production_rate": -1}, {"daughter_expression": 61},
                         {"growth": 1}, {"doubling_time": float("inf")}):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                PopulationConfig(**settings)


if __name__ == "__main__":
    unittest.main()
