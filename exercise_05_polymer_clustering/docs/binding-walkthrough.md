# Follow one binding event

Version 0.2.0 uses the model's three reversible reactions at fixed positions.
Start with `python -m polymer_clustering --mode binding`. This uses the published
bacterium-signal rates with no polymers. Then use the P3 example in the README.

## 1. Separate where an entity is from whether it is bound

`BindingState.from_initial` copies the initial arrays. The original state is
preserved for comparison. It adds bacterial signal occupancy and polymer contact
records, all initially zero or unoccupied. The count of polymer entities does not
change when they acquire attachments.

Suppose a bacterial site contains five free signals and has 15 free receptors.
Its forward propensity is `0.001 * 15 * 5 = 0.075` per arbitrary time unit. This
is a rate, not a probability. Other sites and reaction types have their own rates.

## 2. Construct channels from the current state

`reaction_channels` calculates currently possible events. A bacterial signal
channel aggregates all equivalent receptors at one site. A polymer signal
channel aggregates its five equivalent signal sites. Bacterial contact channels
retain polymer identity so a two-cell bridge can later be broken correctly.

The forward propensity is `k_on * free_sites * available_partners`; the reverse
propensity is `k_off * bound_attachments`. A polymer with two free bacterial ends
has twice the first-attachment propensity of a polymer with one eligible free end,
for otherwise identical partners. Binding capacities are finite.

## 3. Draw a waiting time and select an event

`simulate_binding` adds all propensities to get A. It draws an exponential wait
with mean `1/A`, then selects a channel with probability `channel_rate/A`.
`rng.exponential` performs the waiting-time transformation directly.

If the next event lies beyond the requested duration, it is not applied. If A is
zero, the world remains unchanged to the end of the requested experiment.

## 4. Apply the selected change

For one bacterium-signal binding event:

| Quantity at the site | Before | After |
| --- | --- | --- |
| Free signals | 5 | 4 |
| Bacterium-bound signals | 0 | 1 |
| Free bacterial signal sites | 15 | 14 |
| Total signals in these pools | 5 | 5 |

The reverse propensity now includes `0.002 * 1`. If that unbinding channel is
selected later, one signal returns to the free pool at the same site.

Polymer-signal binding transfers a signal to that polymer's signal count. Polymer-
bacterium binding changes both its contact record and the bacterial site count;
it does not transfer a signal or consume a polymer entity.

## 5. Check the invariants and repeat

`BindingState.validate` checks conservation, capacities, and matching contact
records after every event. A violation raises an error rather than continuing
with an invalid state. The CLI records failure when its output folder is writable.

Observations are taken at user-selected times without altering event times.
`binding_events.csv` contains the actual events; `binding_trajectory.csv` contains
sampled summaries and can miss short-lived changes between observation times.

These outputs describe a binding subsystem. A newly formed bridge changes
connectivity, but does not yet move bacteria or change diffusion. Bacterium-bound
signals are receptor occupancy, not expression or luminescence.
