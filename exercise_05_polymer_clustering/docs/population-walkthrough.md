# Following a division and a production event

Version 0.5.0 adds these processes in `--mode full`. Earlier modes remain useful
for checking smaller pieces. Read POP-01 through POP-06 in the model notes beside
the supplement's S19/S23: some source details require explicit choices.

## A division is an appointment

A triangular distribution describes cell lifetimes in this implementation.
With default mean 15000 and full width 20000 a.u., possible lifetimes run from
5000 to 25000, with the highest density at 15000. We assume the triangle is
symmetric and every initial cell is newborn; those details are not specified
unambiguously by the supplement.

Suppose a cell is assigned a lifetime of 12000 a.u. It has an appointment at
12000, rather than an equal chance to divide during every event. Diffusion and
binding can occur many times before then. If the next stochastic event would
fall after 12000, the solver processes division first. It can then draw a new
exponential wait because the competing constant-rate processes are memoryless.
It does not redraw the scheduled lifetime after every signal hop.

`DivisionSchedule` stores appointments in a heap: the earliest is always next.
Lifetimes and daughter-placement choices have separate seeded random streams.
Observing a map never draws a new lifetime or advances the physical system.

## What changes when a daughter is placed

`nearest_empty_anchor` examines all valid 2 x 2 footprints and finds the nearest
to the parent, using squared Euclidean anchor distance. If several are equally
near, it picks one uniformly. Other bacteria block footprints; small chemicals
do not. A blocked appointment is logged and retried after the configured delay.

On a successful birth:

1. The parent stays in place and retains bound signals, polymers, and expression.
2. The daughter receives a fresh ID and occupies the chosen empty footprint.
3. The daughter has no bound signals or polymers. Free chemicals already at its
   sites remain free and can bind later.
4. Its response clock begins at the birth time, with zero expression by default.
   The initial-cell expression and daughter-expression controls are independent.
5. Both cells receive new division lifetimes, unless the population cap has
   been reached. At the cap, pending division appointments are cleared.
6. The transport/chemical cache is rebuilt because occupancy, receptors,
   obstruction, and singleton-cluster membership have changed.

No bound molecule is copied or divided between cells. That retention rule is
explicit in the source; expression inheritance is not. Our newborn-expression
choice is therefore recorded separately, and may require a sensitivity study.

## BB170 production adds material

At the default BB170 rate, one cell produces on average 0.005 signals per a.u.
Ten cells therefore have a combined propensity of 0.05, giving an average wait
of 20 a.u. when that population stays fixed and production is considered alone.
Individual waits vary; production is a competing Poisson event, not a signal
added at every recorded frame.

When production wins the event competition, `produce_signal` chooses one of the
cells uniformly and one of its four occupied sites uniformly. It adds one free
signal there, then refreshes local binding and diffusion rates. If a division
changes the population, the total production propensity changes immediately.
MM32 produces zero by default; an explicit production-rate override is a model
experiment and is preserved in the settings.

For BB170, starting at zero signals, the invariant is

`free + bacterium-bound + polymer-bound = produced`.

For MM32 the initial 10000 signals are included on the right. Binding moves a
signal between pools; production adds one to both the pool balance and its
expected total. The program checks that balance during and after the run.

## Reading the demonstration

Run `examples/bb170_accelerated_growth_demo.json` to see births quickly. It uses
mean lifetime 40 a.u. and width 20, rather than 15000 and 20000. Plots explicitly
label the changed division timing. It retains the default 0.005 production rate,
so its short duration is not expected to generate a strong quorum response.

In the validated seed-2026 example, 25 cells become 80 by t=150 a.u. There are
55 successful births and 44 new signals. The final pools are 20 free, 20
bacterium-bound, and 4 polymer-bound: 20 + 20 + 4 = 44. The 100 polymers remain
conserved. These are outputs of this demonstration, not published measurements.

Final mean expression is about 0.186 and mean light about 120.69 per cell.
BB170's printed light curve already gives about 114.89 per cell at zero
expression. Much of the population's total light therefore reflects having more
cells, rather than strong expression. Compare the per-cell mean with the sum.
Daughter cells appear near parents even without polymer bridges: proximity alone
is not evidence of polymer-driven clustering. Inspect actual bridge counts and
compare the same model with zero polymers before making a causal claim.

`lineage.csv` connects daughters to parents; `population_trajectory.csv` tracks
births, production, and balance. Frame archives pad unborn cells explicitly;
`present` or `bacteria_counts` distinguishes absence from a real zero response.
GIF playback speed is illustrative; titles give model time in arbitrary units.

## What the assignment still needs

The next scientific comparison is MM32 with and without polymers, using matching
initialisation and observation settings. The supplement's kinetic omissions,
polymer-count/concentration ambiguity, time-axis inconsistency, and light
normalisation remain unresolved; document their effects before claiming agreement.
Then average independent runs with recorded seeds, compare timing, shape and
magnitude against the source figures, and extend to affinity and BB170 cases.
A passing test suite or plausible animation does not replace those comparisons.
