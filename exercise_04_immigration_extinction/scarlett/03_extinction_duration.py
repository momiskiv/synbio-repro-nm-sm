"""Appendix 4: extinction time for q > p using a stable integral representation.
Reconstructed code, not an exact copy of the authors' MATLAB implementation.
"""
from pathlib import Path
import numpy as np
from scipy.special import exp1
from scipy.integrate import quad
import matplotlib.pyplot as plt

OUT = Path(__file__).parent / 'figures'
OUT.mkdir(exist_ok=True)
p, q, dt = 0.4998, 0.5002, 1/24
alpha = q/p
# Extinction-time increments can be expressed as a convergent positive series:
# tau_i = dt/p * sum_{k=0}^infinity alpha^(-k-1)/(i+k)
# Evaluate with quadrature using integral 1/(i+k)=int_0^1 t^(i+k-1)dt.
def increment(i):
    return dt/p * quad(lambda t: t**(i-1) / (alpha-t), 0, 1,
                       epsabs=1e-9, epsrel=1e-9)[0]

sizes = np.arange(1, 1001)
increments = np.array([increment(int(i)) for i in sizes])
durations = np.cumsum(increments)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(sizes, durations, color='black')
ax.set(xlabel='Initial population i', ylabel='Expected extinction time (days)',
       title='Extinction duration (q > p)')
fig.tight_layout()
fig.savefig(OUT / 'extinction_duration.png', dpi=200)
plt.close(fig)
print('Duration at i=1000:', durations[-1])
