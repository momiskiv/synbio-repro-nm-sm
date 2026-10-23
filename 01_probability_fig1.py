"""Reconstructed Figure 1a-c probability curves; not the original executed script."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

OUT = Path(__file__).parent / 'figures'
OUT.mkdir(exist_ok=True)

def establishment_probability(i, N, ratio):
    """Gambler's ruin: chance of hitting N before zero, starting at i."""
    ratio = np.asarray(ratio, dtype=float)
    if np.any(ratio <= 0):
        raise ValueError('q/p must be positive')
    log_r = np.log(ratio)
    with np.errstate(over='ignore', under='ignore', invalid='ignore', divide='ignore'):
        result = np.expm1(i * log_r) / np.expm1(N * log_r)
    result = np.where(np.isclose(log_r, 0, atol=1e-14), i / N, result)
    return result

ratios = np.linspace(0.90, 1.01, 900)
prob = establishment_probability(1, 10_000, ratios)
for suffix, logscale in [('1a', True), ('1b', False)]:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(ratios, prob, color='black')
    if logscale:
        ax.set_yscale('log')
        ax.set_ylim(1e-45, 1)
    ax.set(xlabel='Ratio of deaths to births (q/p)', ylabel='Probability of reaching N = 10,000',
           title=f'Reconstructed Figure {suffix}')
    fig.tight_layout()
    fig.savefig(OUT / f'figure_{suffix}.png', dpi=200)
    plt.close(fig)

N_values = np.unique(np.geomspace(10, 100_000, 400).astype(int))
fig, ax = plt.subplots(figsize=(7, 4.5))
for ratio in (1.000100, 1.000075, 1.000050, 1.000025):
    ax.plot(N_values, [establishment_probability(1, int(n), ratio) for n in N_values],
            label=f'q/p = {ratio:.6f}')
ax.set(xlabel='Target population N', ylabel='Establishment probability',
       title='Reconstructed Figure 1c')
ax.legend()
fig.tight_layout()
fig.savefig(OUT / 'figure_1c.png', dpi=200)
plt.close(fig)
print('Saved reconstructed Figures 1a, 1b and 1c.')
