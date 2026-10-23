# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Font availability, explicit files, all text artists, and scoped plot settings."""

import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
from matplotlib import font_manager
from matplotlib.text import Text
import numpy as np

from polymer_clustering.cli import main
from polymer_clustering.plotting import PlotConfig, plot_style, resolve_font, styled_plot


class PlottingTests(unittest.TestCase):
    def test_missing_requested_font_has_explicit_fallback_or_strict_failure(self) -> None:
        name = 'PolymerTestFontThatDoesNotExist'
        fallback = resolve_font(PlotConfig(font_family=name))
        self.assertTrue(fallback.fallback_used)
        self.assertIn(fallback.family, ('Liberation Sans', 'DejaVu Sans'))
        self.assertEqual(fallback.metadata()['requested_font_family'], name)
        self.assertEqual(fallback.metadata()['resolved_font_family'], fallback.family)
        with self.assertRaises(ValueError):
            resolve_font(PlotConfig(font_family=name, strict_font=True))

    def test_explicit_font_file_and_family_mismatch(self) -> None:
        path = font_manager.findfont(font_manager.FontProperties(family=['DejaVu Sans']),
                                     fallback_to_default=False)
        font = resolve_font(PlotConfig(font_family='DejaVu Sans', font_path=path, strict_font=True))
        self.assertFalse(font.fallback_used)
        self.assertEqual(font.registered_path, Path(path).resolve())
        self.assertEqual(len(font.metadata()['font_file_sha256']), 64)
        with self.assertRaises(ValueError):
            resolve_font(PlotConfig(font_family='Arial', font_path=path))

    def test_all_text_artists_use_selected_font_and_defaults_are_restored(self) -> None:
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        font = resolve_font(PlotConfig(font_family='DejaVu Sans Mono', strict_font=True))
        before = list(matplotlib.rcParams['font.family'])
        with plot_style(font):
            fig, ax = plt.subplots()
            try:
                image = ax.imshow(np.zeros((3, 3)))
                fig.colorbar(image, ax=ax, label='Signals')
                ax.plot([0, 1], [0, 1], label='Example')
                ax.legend()
                ax.set(title='Test title', xlabel='Time', ylabel='Count')
                fig.suptitle('Figure title')
                fig.canvas.draw()
                texts = fig.findobj(match=Text)
                self.assertTrue(texts)
                self.assertTrue(all(t.get_fontfamily() == [font.family] for t in texts))
            finally:
                plt.close(fig)
        self.assertEqual(matplotlib.rcParams['font.family'], before)
        with self.assertRaises(RuntimeError):
            with plot_style(font):
                raise RuntimeError('failed plot')
        self.assertEqual(matplotlib.rcParams['font.family'], before)

    def test_direct_plot_call_uses_policy_without_cli(self) -> None:
        font = resolve_font(PlotConfig(font_family='DejaVu Sans Mono', strict_font=True))
        @styled_plot
        def observe_family():
            return list(matplotlib.rcParams['font.family'])
        with patch('polymer_clustering.plotting.resolve_font', return_value=font):
            self.assertEqual(observe_family(), [font.family])

    def test_invalid_font_requests_fail_before_creating_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)/'bad'
            for flags in (['--font-family', 'PolymerTestFontThatDoesNotExist', '--strict-font'],
                          ['--font-path', str(Path(temporary)/'missing.ttf')]):
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    main([*flags, '--output', str(output)])
                self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
