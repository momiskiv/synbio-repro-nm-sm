# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Check elapsed/CPU distinctions and accounting when a timed phase fails."""

import unittest
from unittest.mock import patch

from polymer_clustering.timing import RunTimer, format_elapsed


class TimingTests(unittest.TestCase):
    def test_independent_wall_cpu_and_phase_measurements(self) -> None:
        with patch('polymer_clustering.timing.perf_counter', side_effect=[10, 11, 14, 20]), \
             patch('polymer_clustering.timing.process_time', side_effect=[2, 7]):
            timer = RunTimer()
            with timer.phase('simulation'):
                pass
            result = timer.snapshot(events=6)
        self.assertEqual(result['total_wall_seconds'], 10)
        self.assertEqual(result['cpu_seconds'], 5)
        self.assertEqual(result['simulation_wall_seconds'], 3)
        self.assertEqual(result['other_wall_seconds'], 7)
        self.assertEqual(result['simulation_events_per_second'], 2)

    def test_failure_records_partial_phase_and_later_phases_accumulate(self) -> None:
        with patch('polymer_clustering.timing.perf_counter', side_effect=[0, 1, 3, 4, 7, 8]), \
             patch('polymer_clustering.timing.process_time', side_effect=[0, 2]):
            timer = RunTimer()
            with self.assertRaises(RuntimeError):
                with timer.phase('output'):
                    raise RuntimeError('rendering failed')
            with timer.phase('output'):
                pass
            result = timer.snapshot()
        self.assertEqual(result['output_wall_seconds'], 5)
        self.assertEqual(result['other_wall_seconds'], 3)
        self.assertIsNone(result['simulation_events_per_second'])
        self.assertIsNone(timer.active_phase)

    def test_long_duration_format(self) -> None:
        self.assertEqual(format_elapsed(62.5), '62.500 s (00:01:02.500)')
        self.assertEqual(format_elapsed(3600), '3600.000 s (01:00:00.000)')


if __name__ == '__main__':
    unittest.main()
