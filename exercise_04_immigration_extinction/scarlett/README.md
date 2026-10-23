# Curtis et al. (2026) — Question 4 Python reproduction

Python scripts reconstructed from our documented analysis of the gambler's ruin microbial immigration/extinction model. **These are reconstructed scripts, not the exact files originally executed during the assignment.** Consequently, plots and numerical details may differ from the figures already submitted; rerun and inspect before claiming exact matches.

## Source
Curtis et al. (2026), *The ISME Journal*, "Probability, parameters, and duration of immigration and extinction in microbial communities".
- Paper: https://academic.oup.com/ismej/article/20/1/wrag195/8740005
- Authors' MATLAB repository: https://github.com/beadyallen/PPD_Immigration_Microbe
- Mathematical reference: appendices 1, 3 and 4 of supplementary material.

## Run
Install Python and run from this folder:

```bash
python -m pip install -r requirements.txt
python 01_probability_fig1.py
python 02_immigration_duration.py
python 03_extinction_duration.py
python 04_supplement_A3_noise.py
python 05_equal_birth_death_validation.py
```

Images are saved in `figures/`. No input data are required.

## What each file does
- `01_probability_fig1.py`: gambler's ruin establishment probabilities, approximating published Figure 1a–c.
- `02_immigration_duration.py`: Appendix 3 recurrence vs direct finite-sum calculation.
- `03_extinction_duration.py`: Appendix 4 extinction-duration calculation (q>p).
- `04_supplement_A3_noise.py`: illustrative random death-probability perturbations (seed 42).
- `05_equal_birth_death_validation.py`: Appendix 4 equal birth/death carrying-capacity boundary comparison.

## Limitations and transparency
- Figure 1c and noise trajectories are qualitative reconstructions, not independently validated point-for-point against authors' MATLAB results.
- The published Figure A3 caption states a +/-0.0002 perturbation; the authors' MATLAB noise setting may differ. This reconstruction samples one q value per trajectory rather than necessarily matching the author's per-step noise procedure.
- Immigration duration follows the paper's recurrence for the unconditioned hitting-time formulation; care is required interpreting time-to-target when extinction can occur first.
- The equal-probability boundary comparison is a validation experiment; it is not a direct copy of a published figure.
- Both group members worked on the same computer; independent rerunning on separate machines has **not** been completed.
- Python code was reconstructed with ChatGPT assistance from the documented equations and parameters; review the outputs before submission.

## Previously reported environment (verify on your own computer)
- Python 3.15.0
- NumPy 2.5.4
- Matplotlib 3.11.2
- SciPy 1.18.1

## AI prompts
See `Q4_AI_Prompts.md` for a paraphrased record of Python-related AI prompts and disclosure of AI assistance.
