"""Illustrative reconstruction of Figure A3 noise behaviour.
The published caption and MATLAB code use different noise specifications;
this script follows the caption's +/-0.0002 range, not exact MATLAB RNG.
"""
from pathlib import Path
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

OUT = Path(__file__).parent / 'figures'
OUT.mkdir(exist_ok=True)
DT = 1/24
Q_BASE = 0.5002
SEED = 42
N_RUNS = 20
N_MAX = 1000

def extinction_curve(q):
    p = 1-q
    if q <= p:
        return None  # unbounded extinction expectation not finite
    alpha = q/p
    increments = [DT/p * quad(lambda t: t**(i-1)/(alpha-t), 0, 1)[0]
                  for i in range(1, N_MAX+1)]
    return np.cumsum(increments)

rng = np.random.default_rng(SEED)
sizes = np.arange(1, N_MAX+1)
fig, ax = plt.subplots(figsize=(7, 4.5))
for _ in range(N_RUNS):
    q = Q_BASE + rng.uniform(-0.0002, 0.0002)
    curve = extinction_curve(q)
    if curve is not None:
        ax.plot(sizes, curve, color='red', alpha=0.28, lw=0.8)
ax.plot(sizes, extinction_curve(Q_BASE), color='black', lw=2, label='Baseline q=0.5002')
ax.set(xlabel='Initial population i', ylabel='Expected extinction time (days)',
       title='Supplementary Figure A3: illustrative noise reconstruction')
ax.legend()
fig.tight_layout()
fig.savefig(OUT / 'supplement_A3_noise.png', dpi=200)
plt.close(fig)
