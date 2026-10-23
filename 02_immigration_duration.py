"""Appendix 3: recurrence and direct summation of expected time to reach n.
Reconstructed script. These formulas do not condition on eventual success.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

OUT = Path(__file__).parent / 'figures'
OUT.mkdir(exist_ok=True)
p, q, dt = 0.4998, 0.5002, 1/24
N_MAX = 1000

def recurrence(n):
    step_time = 0.0
    total = 0.0
    for i in range(1, n):
        step_time = dt / (i * p) + (q / p) * step_time
        total += step_time
    return total

def direct_sum(n):
    # Equation 3.5 written as finite sums; slower but independent of recurrence.
    alpha = q / p
    return sum(dt / p * sum(alpha ** (i - j) / j for j in range(1, i + 1))
               for i in range(1, n))

sizes = np.arange(1, N_MAX + 1)
rec = np.array([recurrence(int(n)) for n in sizes])
# Efficient direct summation of the same finite-series formula.
alpha = q / p
step_direct = np.zeros(N_MAX)
for i in range(1, N_MAX):
    # finite sum in Eq 3.5 (computed without recursive step_time)
    step_direct[i] = dt / p * sum(alpha ** (i-j) / j for j in range(1, i+1))
ana = np.cumsum(step_direct)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(sizes, rec, label='Recurrence', lw=2)
ax.plot(sizes, ana, '--', label='Direct finite sum', lw=1.5)
ax.set(xlabel='Target population n', ylabel='Expected duration (days)', title='Immigration duration')
ax.legend()
fig.tight_layout()
fig.savefig(OUT / 'immigration_duration.png', dpi=200)
plt.close(fig)
print('Max absolute difference:', np.max(np.abs(rec - ana)))
