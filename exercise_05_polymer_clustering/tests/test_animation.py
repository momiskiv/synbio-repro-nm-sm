# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Regression check for colourbar geometry changing during animation recording."""

import unittest

import matplotlib
import numpy as np

from polymer_clustering.output import freeze_animation_layout


class AnimationLayoutTests(unittest.TestCase):
    def test_colourbar_bounds_and_numeric_limits_stay_fixed(self) -> None:
        matplotlib.use("Agg")
        from matplotlib import pyplot as plt

        figure, axes = plt.subplots(1, 2, figsize=(11, 5), layout="constrained")
        try:
            images = []
            for ax in axes:
                image = ax.imshow(np.zeros((50, 50)), vmin=0, vmax=16, origin="lower")
                figure.colorbar(image, ax=ax, label="Fixed units", shrink=0.75)
                images.append(image)
            figure.suptitle("Model response | t = 0 a.u.")
            freeze_animation_layout(figure)
            expected = [ax.get_window_extent().bounds for ax in figure.axes]
            for time in (0, 5, 50, 100):
                figure.suptitle(f"Model response | t = {time} a.u.")
                for image in images:
                    image.set_data(np.full((50, 50), time % 17))
                figure.canvas.draw()
                actual = [ax.get_window_extent().bounds for ax in figure.axes]
                np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-7)
                self.assertTrue(all(image.get_clim() == (0, 16) for image in images))
        finally:
            plt.close(figure)


if __name__ == "__main__":
    unittest.main()
