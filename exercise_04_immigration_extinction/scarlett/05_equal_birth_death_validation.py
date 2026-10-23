"""Appendix 4, p=q=0.5: compare carrying-capacity boundary implementations.
Reconstructed validation, not original executed script.
"""
from pathlib import Path
import numpy as np
from scipy.special import digamma
import matplotlib.pyplot as plt

OUT = Path(__file__).parent / 'figures'
OUT.mkdir(exist_ok=True)
DT = 1/24
SIZES = np.arange(1, 1001)
EULER_GAMMA = np.euler_gamma

def harmonic(n):
    return digamma(n+1) + EULER_GAMMA

def analytic(i, capacity):
    # Equation 4.17
    return 2*DT*i*(harmonic(capacity-1)-harmonic(i-1)) + i/capacity + 2*DT*(i-1)

def recurrence(capacity, boundary_capacity):
    # Difference in boundary condition: authors' code reportedly uses 1000.
    tau = 1/boundary_capacity + 2*DT*harmonic(capacity-1)
    total = 0.0
    values = []
    for i in SIZES:
        if i > 1:
            tau -= 2*DT/(i-1)
        total += tau
        values.append(total)
    return np.array(values)

fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, capacity in zip(axs, (10_000, 10**16)):
    exact = np.array([analytic(int(i), capacity) for i in SIZES])
    original_boundary = recurrence(capacity, 1000)
    corrected = recurrence(capacity, capacity)
    ax.plot(SIZES, exact, label='Analytical', lw=2)
    ax.plot(SIZES, original_boundary, '--', label='Recurrence: 1/1000 boundary')
    ax.plot(SIZES, corrected, ':', label='Recurrence: 1/capacity boundary')
    ax.set(title=f'Capacity = {capacity:.0e}', xlabel='Initial population i',
           ylabel='Expected extinction time (days)')
    print(f'Capacity {capacity}: corrected maximum error {np.max(np.abs(exact-corrected)):.3e}')
axs[0].legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / 'equal_birth_death_validation.png', dpi=200)
plt.close(fig)
