# From receptor occupancy to predicted light

Use `--mode response`, or the supplied `mm32_p3_response.json` example. This mode
adds a downstream response to the same binding and movement simulation. The
source is supplementary S21-S22; interpretation choices are RESP-01 through
RESP-05 in model-notes.md.

## Three values with different meanings

| Value per cell | Meaning | How it changes |
| --- | --- | --- |
| Sensing S | Signals bound to that bacterium, summed over its four sites | Jumps by one on a BS binding/release event |
| Expression E | Simplified downstream expression state | Changes continuously with elapsed time |
| Light L | Predicted response derived from E | Follows the nonlinear S22 transformation |

Polymer-bound signals do not contribute to S. Free signals sharing a cell's site
also do not contribute until they bind. Sensing is an integer, but expression
and predicted light are continuous model values.

## Why expression lags behind sensing

If S rises to 20 while E is 5, the instantaneous synthesis rate is
`0.05 * (20 - 5) = 0.75`. That is a rate per time unit, not an amount to add after
each stochastic event. As the gap shrinks, synthesis slows.

If S stays at 20 for the next 20 units, the exact solution gives
`E = 20 + (5 - 20)*exp(-0.05*20)`, approximately 14.48. It has moved towards 20
but has not jumped straight there. The synthesis time constant is 1/0.05=20 a.u.

If S falls below E, the specified expression variable decays with positive rate
magnitude 0.01, using the printed leading minus sign. At fixed nonzero S, decay
stops when E reaches S, because the falling condition no longer holds. At S=0,
decay approaches zero exponentially. This branch handling is documented rather
than silently applying a different relaxation equation.

## How the code handles elapsed time

`ResponseTracker` stores each cell's sensing level, expression value, and the
time at which its sensing last changed. Before a new BS event changes sensing,
it evaluates expression up to that event under the old sensing level.

For a plot or GIF frame, it evaluates the same exact solution at the requested
time without changing the stored integration anchors. Receptor sensing does not
change during unrelated signal hops, so they do not require all cells to be
updated. This is both more efficient and independent of observation frequency.
Movement changes a cell's position while retaining its response identity.

## Light is a smooth response, not a switch

The printed transformation is `L = 2*Scale / (1 + exp(0.25*(BSc-E)))`.
At E=BSc, the exponential is 1 and L=Scale. For MM32, BSc=38 and Scale=2500.
An E below 38 can still produce nonzero light; it is not an on/off threshold.

The formula also gives a small nonzero value at E=0. We retain that baseline.
The theoretical upper asymptote is 5000, although the allowed expression range
does not necessarily reach it. These numbers are model arbitrary units, not
measured photons. BB170 uses a different BSc. Response diagnostic mode excludes
production; full mode adds it, as described in the population walkthrough.

Each cell's light is calculated before population aggregation. For cells with
different expression levels, applying the nonlinear function to their average
would give a different answer and lose population variation.

## Read the plots and animations together

`response.png` shows mean sensing and expression on the left, and mean predicted
light on the right. `response.gif` pairs free-signal distribution with each cell's
predicted light. A locally bright signal patch need not immediately create a
bright cell: free signal must bind, and expression then takes time to respond.

The light colourbar normally spans this run's observed values, with a minimum
upper display bound of 1. A bright colour can therefore represent a dim absolute
response. Use the same `--light-map-max` for visual comparisons between runs.
The per-cell and population CSVs preserve exact unscaled model values.

Version 0.5.0 adds division and production in full mode; see
[the population walkthrough](population-walkthrough.md).
Do not interpret a short illustrative run as the full published time course.
