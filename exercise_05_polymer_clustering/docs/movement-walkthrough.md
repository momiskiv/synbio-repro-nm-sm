# Follow a movement event

The new mode is `--mode movement`. It adds transport events to binding in the
same Gillespie process. It still has no growth, signal production, or light.

## A free signal hops

Suppose a site contains four free signals, no polymers, and no bacterial
obstruction. Its propensity to send one signal to a particular neighbouring
site is `4 * D_S * f(0)`. D_S is 0.05 and the exact printed mobility function gives
f(0) approximately 0.999695. Each of the four legal directions has its own rate.

When one hop is selected, the source loses one free signal and the destination
gains one. The total is unchanged. The source and destination binding rates are
then updated: both the local signal counts and opportunities for receptor
binding have changed. Bound signals cannot independently make these hops.

At a closed edge, directions outside the grid have zero propensity. Directions
touching a bacterial footprint are multiplied by the chosen obstruction factor.
That factor is configurable because the source does not specify an exact value.

## A polymer moves with its captured signals

Each polymer has an identity, coordinate, signal occupancy, and up to two
bacterial contacts. Only a polymer with no bacterial contact diffuses on its own.
Its bound signals remain attached during a hop. It can subsequently release them
at its new location. This is distinct from moving free signals.

Local polymer counts affect signal mobility and the accessibility of polymer
signal-binding sites. Both must be updated after a polymer hop.

## Bacteria move together when bridged

Connected components are built from individual polymer bridges. Merely adjacent
bacteria remain separate components. For a two-cell component, the per-direction
movement rate is `0.0005 / 2 = 0.00025`, if every member can move in that direction.

A selected translation changes all member coordinates simultaneously. Bacterial
receptor occupancy, tethered polymers, and both bridge contacts move with them.
Another bacterium or an outer edge can block the whole component. Releasing a
bridge can split a component, so components must be recalculated after PB events.

Our documented displacement rule also shifts free chemicals in the old/new
footprint region along the travel direction, conserving them at closed edges.
This is an interpretation of the supplement, not a detailed fluid-flow model.

## Why the solver stores cached rates

With 10000 free signals and a per-direction rate near 0.05, there can initially
be roughly 2000 signal hops per time unit before boundary/obstruction reductions.
Rebuilding every channel after each hop would be costly.

The solver first selects an event family, then a site/entity within that family,
then a direction or chemical channel. Selecting a family with rate A_family out
of total A, and a channel with rate a out of A_family, gives `(A_family/A) *
(a/A_family) = a/A`, the ordinary Gillespie probability. The waiting time uses
the total rate A. This is still event-based SSA, not a fixed time-step simulation.

Most signal hops affect only two sites. Bacterial motion and PB contact changes
trigger a full geometry refresh. Tests compare cached rates to full recalculation
after mixed events and compare unobstructed diffusion to its analytical law.

## Read the new outputs

- `binding.png` shows the three signal pools and attachment counts over time.
- `movement.png` compares initial and final free-signal distributions.
- `movement.gif` shows spatial observations, with simulation time in each frame.
- `movement_positions.csv` and `final_bacteria.csv` record cell coordinates.
- `movement_frames.npz` stores the plotted arrays and positions.

Frame rate affects playback, not the simulation. A 100-a.u. run may show only a
few bacterial translations because the source bacterial rate is much lower than
the signal rate. A run need not contain a bridge; bridge motion is separately
tested with controlled two-cell states. These outputs do not yet show light.
