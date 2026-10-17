#Module 1: Synthetic Biology and Bioinformatics - EngBio4Env
#Working on Curtis et al. (2026) - Probabilities of survival within microbial communities

#Author: N.Musto
#Last update: 17-10-2026

import argparse #parameters
import math #for log scaling
import matplotlib.pyplot as plt #graphs!
from matplotlib import font_manager
from matplotlib.ticker import (
    MultipleLocator, FormatStrFormatter, 
    MaxNLocator, AutoMinorLocator) #formatting graphs
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable #colour bar for Fig1c
import csv
from scipy.integrate import quad #evaluation of definite integral and returns estimate numerical error
import platform
from importlib.metadata import version #versions; for later writeup and for csv
import numpy as np #local seeded generator
from datetime import datetime #self-explanatory
from pathlib import Path #file system paths
import tempfile #folder creation


def log10_target_prob(i, N, ratio):
    if N < 1 or not 0 <= i <= N:
        raise ValueError("Require N >= 1 and 0 <= i <= N.")

    if ratio < 0 or not math.isfinite(ratio):
        raise ValueError("The ratio must be finite and non-negative.")

    # Boundary cases
    if i == 0:
        return -math.inf

    if i == N or ratio == 0:
        return 0.0

    # Equal birth and death probabilities
    if ratio == 1.0:
        return math.log10(i) - math.log10(N)

    log_r = math.log(ratio)

    if ratio < 1:
        log_prob = (
            math.log(-math.expm1(i * log_r))
            - math.log(-math.expm1(N * log_r)))
    else:
        log_prob = (
            (i - N) * log_r
            + math.log(-math.expm1(-i * log_r))
            - math.log(-math.expm1(-N * log_r)))

    return log_prob / math.log(10)
#Added because plotting large exponentials could be an issue later with larger numbers



def target_prob(i, N, ratio):
    return 10**log10_target_prob(i, N, ratio)
    #As seen on the supplementary material of the paper; equation 1
    

def save_probability_csv(filename, figure, initial, targets, ratios):
    with open(filename, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)

        writer.writerow([
            "figure", "initial", "target", "ratio_q_over_p",
            "p", "q", "probability", "log10_probability",
            "log10_target", "probability_underflow"])

        for current_ratio in ratios:
            birth = 1 / (1 + current_ratio)
            death = 1 - birth

            for current_target in targets:
                log_probability = log10_target_prob(
                    initial, current_target, current_ratio)
                probability = 10**log_probability

                underflow = (
                    probability == 0
                    and math.isfinite(log_probability))

                writer.writerow([
                    figure, initial, current_target, current_ratio,
                    birth, death, probability, log_probability,
                    math.log10(current_target), underflow])

    print(f"Saved probability data: {filename}")

def extinction_times(max_initial, ratio, delta_t): #originally inside a MATLAB script #check out the gambler's ruin on their repo
    """Return T_0 to T_max_initial, in the units of delta_t."""

    if max_initial < 0:
        raise ValueError("Starting abundance must be non-negative.")

    if not math.isfinite(ratio) or ratio <= 1:
        raise ValueError(
            "This extinction calculation requires q/p > 1.")

    if not math.isfinite(delta_t) or delta_t <= 0:
        raise ValueError("The time scale must be finite and positive.")

    # T_0 = 0: extinction has already occurred.
    times = [0.0] #accumulated extinction times

    if max_initial == 0:
        return times

    birth = 1 / (1 + ratio)

    # tau_1 = T_1
    increment = -(delta_t / birth) * math.log1p(-1 / ratio) #t_i
    times.append(increment)

    for current_i in range(1, max_initial):
        # tau_(i+1) = ratio * tau_i - delta_t / (p * i)
        increment = (
            ratio * increment
            - delta_t / (birth * current_i))

        if not math.isfinite(increment) or increment <= 0:
            raise ValueError(
                "Unstable forward recurrence for these settings.")

        # T_(i+1) = T_i + tau_(i+1)
        times.append(times[-1] + increment)

    return times

def extinction_time_analytic(initial, ratio, delta_t):
    """Return mean extinction time and estimated quadrature error."""

    if initial < 0:
        raise ValueError("Starting abundance must be non-negative.")

    if not math.isfinite(ratio) or ratio <= 1:
        raise ValueError(
            "This extinction calculation requires q/p > 1.")

    if not math.isfinite(delta_t) or delta_t <= 0:
        raise ValueError("The time scale must be finite and positive.")

    if initial == 0:
        return 0.0, 0.0

    birth = 1 / (1 + ratio)

    def integrand(x): #integration
        if x == 0:
            geometric_sum = 1.0
        elif x == 1:
            geometric_sum = float(initial)
        else:
            geometric_sum = ( 
                -math.expm1(initial * math.log(x)) / (1 - x))

        return geometric_sum / (ratio - x)

    integral, error = quad(
        integrand,
        0.0,
        1.0,
        epsabs=1e-10,
        epsrel=1e-10,
        limit=200)

    scale = delta_t / birth

    return scale * integral, scale * error

def growth_times(max_target, ratio, delta_t): #Figure A1
    """Return {n: M_1,n} for targets 1 to max_target, using m_0 = 0."""

    if max_target < 1:
        raise ValueError("Target abundance must be at least 1.")

    if not math.isfinite(ratio) or ratio < 0:
        raise ValueError(
            "The ratio must be finite and non-negative.")

    if not math.isfinite(delta_t) or delta_t <= 0:
        raise ValueError(
            "The time scale must be finite and positive.")

    birth = 1 / (1 + ratio)

    # Starting at one cell already meets a target of one.
    times = {1: 0.0} #target abundance as keys

    # The appendix's starting condition: m_0 = 0.
    increment = 0.0 #m_i

    for current_i in range(1, max_target):
        increment = (
            delta_t / (birth * current_i)
            + ratio * increment)

        total = times[current_i] + increment

        if not math.isfinite(total):
            raise ValueError(
                "The growth-time recurrence overflowed for these settings.")

        times[current_i + 1] = total

    return times

def growth_time_analytic(target, ratio, delta_t):
    """Evaluate the A1 finite-sum solution with m_0 = 0."""

    if target < 1:
        raise ValueError("Target abundance must be at least 1.")

    if not math.isfinite(ratio) or ratio < 0:
        raise ValueError(
            "The ratio must be finite and non-negative.")

    if not math.isfinite(delta_t) or delta_t <= 0:
        raise ValueError(
            "The time scale must be finite and positive.")

    if target == 1:
        return 0.0

    birth = 1 / (1 + ratio)

    try:
        if ratio == 0:
            # Pure birth: no deaths.
            total = math.fsum( #adds terms with improved floating point accuracy
                1 / k for k in range(1, target))

        elif ratio == 1:
            # Equal birth and death probabilities.
            total = math.fsum(
                (target - k) / k
                for k in range(1, target))

        else:
            log_ratio = math.log(ratio)

            total = math.fsum(
                math.expm1((target - k) * log_ratio) #evaluates the numerator when ratio nears 1
                / ((ratio - 1) * k)
                for k in range(1, target))

        result = (delta_t / birth) * total #mean duration

    except OverflowError:
        raise ValueError(
            "The analytical growth sum overflowed for these settings."
        ) from None

    if not math.isfinite(result):
        raise ValueError(
            "The analytical growth sum overflowed for these settings.")

    return result

def prepare_output_prefix(figure, output_dir="results", output_prefix=None):
    """Create a figure destination and return its shared filename prefix."""

    if output_prefix is not None:
        if not output_prefix.strip():
            raise ValueError("The output prefix must not be empty.")

        prefix = Path(output_prefix).expanduser().resolve()
        prefix.parent.mkdir(parents=True, exist_ok=True)

    else:
        root = Path(output_dir).expanduser().resolve()
        root.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().astimezone().strftime(
            "%Y%m%d_%H%M%S")

        run_dir = Path(tempfile.mkdtemp(
            prefix=f"{timestamp}_figure{figure}_", dir=root))

        prefix = run_dir / f"figure{figure}"

    print(f"\nOutput folder: {prefix.parent}")

    return str(prefix)

#Flags

parser = argparse.ArgumentParser(
    description="Calculate the probability of reaching a target abundance.")

parser.add_argument("--initial", type=int, default=1,
                    help="Starting population (default: 1)")

parser.add_argument("--target", type=int, default=10**4,
                    help="Target abundance (default: 10**4)")

parser.add_argument("--ratio", type=float, default=1.0,
                    help="Death-to-birth ratio q/p (default: 1.0)")

parser.add_argument("--ratios", type=float, nargs="+",
                    help="Multiple death-to-birth ratios to compare in Figure 1c")

parser.add_argument("--figure", choices =["1a", "1b", "1c", "A1", "A2", "A3"],
                    help="Generate a selected figure (available: ALL")

parser.add_argument("--duration", choices=["extinction", "growth"],
    help="Calculate a selected mean duration")

parser.add_argument("--time-scale", type=float, default=1 / 24,
    help="Time constant delta_t in days (default: 1/24)")

parser.add_argument("--method", choices=["recurrence", "analytic", "both"], default="recurrence",
    help="Duration calculation method (default: recurrence)")

parser.add_argument("--max-initial", type=int, default=1000,
    help="Largest starting abundance plotted in A2/A3 (default: 1000)") #selects upper end of the plotted range

parser.add_argument("--max-target", type=int, default=1000,
    help="Largest target abundance plotted in A1 (default: 1000)")

parser.add_argument("--noise-halfwidth", type=float, default=0.0002,
    help="Half-width of uniform noise added to q in A3")

parser.add_argument("--draws", type=int, default=20,
    help="Number of A3 parameter draws (default: 20)")

parser.add_argument("--show-draws", type=int, default=20,
    help="Maximum number of parameter curves displayed")

parser.add_argument("--seed", type=int, default=2026,
    help="Random seed for A3")

parser.add_argument("--summary-only", action="store_true",
    help="For A3, skip the full curve CSV; retain summaries and plots")

parser.add_argument("--output-dir", default="results",
    help="Parent for automatic figure folders (default: ./results)")

parser.add_argument("--output-prefix", default=None,
    help="Exact output path prefix; overrides automatic figure folders")

args = parser.parse_args()

if args.ratios is not None and args.figure != "1c":
    parser.error("--ratios requires --figure 1c.")

if args.duration and args.figure:
    parser.error(
        "Choose either --duration or --figure for this run.")

if args.figure == "A3" and args.method != "recurrence":
    parser.error("A3 currently uses --method recurrence.")

if args.method != "recurrence" and not (
    args.duration in ["extinction", "growth"]
    or args.figure in ["A1", "A2"]):
    parser.error(
        "--method analytic/both requires a duration calculation or A1/A2.")

if args.figure in ["A1", "A2"] and args.method != "both":
    parser.error(
        "Figures A1 and A2 compare both methods: add --method both.")

if args.duration == "growth" or args.figure == "A1":
    if args.initial != 1:
        parser.error(
            "A1 growth calculations require --initial 1.")

#Mathematical Variables

i = args.initial
N = args.target
ratio = args.ratio

if not math.isfinite(ratio) or ratio < 0:
    parser.error(
        "The ratio must be finite and non-negative.")
    
p = 1 / (1+ ratio)
q = 1 - p

#CLI Output

if args.duration == "growth":
    try:
        if args.method in ["recurrence", "both"]:
            growth = growth_times(
                N, ratio, args.time_scale)

        if args.method in ["analytic", "both"]:
            analytic_growth = growth_time_analytic(
                N, ratio, args.time_scale)

    except ValueError as error:
        parser.error(str(error))

    print("\nA1 growth-duration calculation:")
    print("Starting population: 1")
    print(f"Target abundance: {N}")
    print("Lower boundary: holding at 1 (m_0 = 0)")
    print(f"Death-to-birth ratio (q/p): {ratio}")
    print(f"Birth probability: {p}")
    print(f"Death probability: {q}")
    print(f"Time constant delta_t: {args.time_scale} days")
    print(f"Method: {args.method}")

    if args.method in ["recurrence", "both"]:
        print(f"Recurrence mean: {growth[N]:.9f} days")
        print(f"Recurrence mean: {24 * growth[N]:.9f} hours")

    if args.method in ["analytic", "both"]:
        print(f"Analytical mean: {analytic_growth:.9f} days")
        print(f"Analytical mean: {24 * analytic_growth:.9f} hours")

    if args.method == "both":
        difference = abs(growth[N] - analytic_growth)
        print(f"Absolute difference: {difference:.3e} days")

    parser.exit()

if args.duration == "extinction":
    try:
        if args.method in ["recurrence", "both"]:
            times = extinction_times(
                i, ratio, args.time_scale)

        if args.method in ["analytic", "both"]:
            analytic_time, quadrature_error = extinction_time_analytic(
                i, ratio, args.time_scale)

    except ValueError as error:
        parser.error(str(error))

    print("\nMean time to extinction:")
    print(f"Starting population: {i}")
    print("Stopping abundance: 0")
    print("Upper population boundary: none")
    print(f"Death-to-birth ratio (q/p): {ratio}")
    print(f"Birth probability: {p}")
    print(f"Death probability: {q}")
    print(f"Time constant delta_t: {args.time_scale} days")
    print(f"Method: {args.method}")

    if args.method in ["recurrence", "both"]:
        print(f"Recurrence mean: {times[-1]:.9f} days")
        print(f"Recurrence mean: {24 * times[-1]:.9f} hours")

    if args.method in ["analytic", "both"]:
        print(f"Analytical mean: {analytic_time:.9f} days")
        print(f"Analytical mean: {24 * analytic_time:.9f} hours")
        print(
            f"Estimated quadrature error: "
            f"{quadrature_error:.3e} days"
        )

    if args.method == "both":
        difference = abs(times[-1] - analytic_time)
        print(f"Absolute difference: {difference:.3e} days")

    parser.exit()

if args.figure not in ["A1", "A2", "A3"]:
    print(
        "\nSingle-target calculation using "
        "--initial, --target and --ratio:")
    print(f"Starting population: {i}")
    print(f"Target abundance: {N}")
    print(f"Death-to-birth ratio (q/p): {ratio}")
    print(f"Birth probability: {p}")
    print(f"Death probability: {q}")
    print(f"Probability of reaching target: {target_prob(i, N, ratio)}")
    print(f"Log10 probability: {log10_target_prob(i, N, ratio)}")

#Figures (pdf, svg, png)

#Figure 1a&b

font_manager.findfont("Arial", fallback_to_default=False) #I like the look of Arial for the figure

if args.figure:
    plt.rcParams.update({
            "font.family": "Arial",
            "font.size": 9,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.7,

            # Mathematical notation
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "mathtext.sf": "Arial",
            "mathtext.cal": "STIXGeneral",
            "mathtext.fallback": "stix",

            # Export settings
            "pdf.fonttype": 42,
            "svg.fonttype": "path",
        })

# Resolve the destination once, before generating the selected figure.
if args.figure:
    try:
        args.output_prefix = prepare_output_prefix(
            args.figure, args.output_dir, args.output_prefix)

    except (ValueError, OSError) as error:
        parser.error(str(error))

if args.figure in ["1a", "1b"]:
    ratios = []
    plot_values = []

    for step in range(920, 1041):
        current_ratio = step / 1000
        probability = target_prob(i, N, current_ratio)

        if args.figure == "1a":
            plot_value = log10_target_prob(i, N, current_ratio)
        else:
            plot_value = probability

        ratios.append(current_ratio)
        plot_values.append(plot_value)


    fig, ax = plt.subplots(
        figsize=(3.7, 3.5),
        layout="constrained"
    )

    ax.plot(
        ratios, plot_values,
        color="black",
        linewidth=2.0,
        marker="o",
        markersize=1.5,
        solid_capstyle="round",
        zorder=3
    )

    ax.axvline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=0.9,
        zorder=2
    )

    ax.set_xlabel(
        r"$\left(\frac{q}{p}\right)$ Ratio of deaths to births",
        labelpad=7
    )

        # Format powers of ten as exponents; show other targets normally.
    exponent = len(str(N)) - 1

    if N == 10**exponent:
        target_label = rf"10^{{{exponent}}}"
    else:
        target_label = str(N)

    if args.figure == "1a":
        vertical_label = r"$\log_{10}(P_i)$"
    else:
        vertical_label = r"$P_i$"

    ax.set_ylabel(
        vertical_label + "\n"
        + rf"Probability of reaching ${target_label}$ individuals",
        labelpad=8
    )

    # Small margins keep the curve away from the frame.
    ax.set_xlim(0.915, 1.045)
    if args.figure == "1a":
        ax.set_ylim(-180, 8)
        ax.yaxis.set_major_locator(MultipleLocator(50))
        ax.yaxis.set_minor_locator(MultipleLocator(25))
    else:
        upper_limit = max(plot_values) * 1.06

        ax.set_ylim(-0.025 * upper_limit, upper_limit)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_facecolor("#FCFCFC")

    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_minor_locator(MultipleLocator(0.025))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#777777")
        spine.set_linewidth(0.7)

    ax.tick_params(
        which="major", direction="out",
        length=3, width=0.7
    )
    ax.tick_params(which="minor", length=0)

    ax.set_axisbelow(True)
    ax.grid(which="major", color="#AFAFAF", linewidth=0.65)
    ax.grid(which="minor", color="#D6D6D6", linewidth=0.45)

    output_name = args.output_prefix
    save_probability_csv(
        filename=f"{args.output_prefix}.csv",
        figure=args.figure,
        initial=i,
        targets=[N],
        ratios=ratios
    )

    fig.savefig(f"{args.output_prefix}.png", dpi=600, facecolor="white")
    fig.savefig(f"{args.output_prefix}.pdf", facecolor="white")
    fig.savefig(f"{args.output_prefix}.svg", facecolor="white")

    plt.show()

#Figure 1c

if args.figure == "1c":
    # Use multiple ratios when supplied; otherwise keep single-ratio mode.
    if args.ratios is None:
        curve_ratios = [ratio]
    else:
        curve_ratios = sorted(set(args.ratios), reverse=True) #from largest

    print(f"\nFigure 1c curve ratios: {curve_ratios}")

    targets = []
    log_targets = []

    # Targets spaced by 0.1 on the log10 scale.
    for step in range(10, 81):
        current_target = round(10**(step / 10))

        if current_target < i:
            continue

        targets.append(current_target)
        log_targets.append(math.log10(current_target))

    if not targets:
        parser.error("The starting population exceeds the target range.")

    fig, ax = plt.subplots(
        figsize=(3.7, 3.5),
        layout="constrained"
    )

    # The colour scale represents marker fills.
    cmap = LinearSegmentedColormap.from_list(
        "ratio_palette",
        ["#5041AB", "#C146A2"]
    )

    if len(curve_ratios) > 1:
        norm = Normalize(
            vmin=min(curve_ratios),
            vmax=max(curve_ratios)
        )

    for current_ratio in curve_ratios:
        log_probabilities = []

        for current_target in targets:
            log_probabilities.append(
                log10_target_prob(i, current_target, current_ratio)
            )

        if len(curve_ratios) > 1:
            face_colour = cmap(norm(current_ratio))

            # Slightly darker outlines give subtle contrast.
            edge_colour = tuple(
                0.85 * channel for channel in face_colour[:3]
            )
        else:
            face_colour = "#C146A2"
            edge_colour = "#B2188B"

        ax.plot(
            log_targets,
            log_probabilities,
            linestyle="none",
            marker="o",
            markersize=3.2,
            markerfacecolor=face_colour,
            markeredgecolor=edge_colour,
            markeredgewidth=0.9,
            label=rf"$\frac{{q}}{{p}}={current_ratio}$",
            zorder=3
        )

    ax.set_xlabel(r"$(N)$ Population size", labelpad=7)
    ax.set_ylabel(r"$P_i$", labelpad=8)

    ax.set_xlim(0.75, 8.25)
    ax.set_ylim(-52, 2)
    ax.set_facecolor("#FCFCFC")

    ax.xaxis.set_major_locator(MultipleLocator(2.5))
    ax.xaxis.set_minor_locator(MultipleLocator(1.25))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))

    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(5))

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#777777")
        spine.set_linewidth(0.7)

    ax.tick_params(
        which="major", direction="out",
        length=3, width=0.7
    )
    ax.tick_params(which="minor", length=0)

    ax.set_axisbelow(True)
    ax.grid(which="major", color="#AFAFAF", linewidth=0.65)
    ax.grid(which="minor", color="#D6D6D6", linewidth=0.45)

    if len(curve_ratios) > 1:
        # Position and size, expressed as fractions of the plotting area.
        cax = ax.inset_axes([0.05, 0.07, 0.035, 0.30])

        colour_bar = fig.colorbar(
            ScalarMappable(norm=norm, cmap=cmap),
            cax=cax
        )

        colour_bar.set_label(
            r"$\frac{q}{p}$",
            rotation=0,
            labelpad=8,
            fontsize=9
        )

        colour_bar.locator = MaxNLocator(nbins=4)
        colour_bar.formatter = FormatStrFormatter("%.6f")
        colour_bar.update_ticks()

        colour_bar.ax.tick_params(
            labelsize=7,
            length=2,
            width=0.6
        )
        colour_bar.outline.set_linewidth(0.6)

        ax.text(
            0.96, 0.96,
            rf"$i={i}$",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8
        )
    else:
        ax.legend(
            loc="lower left",
            frameon=False,
            fontsize=8,
            title=rf"$i={i}$",
            title_fontsize=8
        )

    save_probability_csv(
        filename="figure1c.csv",
        figure="1c",
        initial=i,
        targets=targets,
        ratios=curve_ratios
    )

    fig.savefig("figure1c.png", dpi=600, facecolor="white")
    fig.savefig("figure1c.pdf", facecolor="white")
    fig.savefig("figure1c.svg", facecolor="white")

    plt.show()

# Figure A1: comparison of reported growth-duration calculations

if args.figure == "A1":
    if args.max_target < 2:
        parser.error(
            "The A1 log plot requires --max-target of at least 2.")

    settings = {
        "initial": 1,
        "ratio_q_over_p": ratio,
        "p": p,
        "q": q,
        "delta_t_days": args.time_scale,
        "lower_boundary": "m_0=0; equivalent holding at 1",
        "upper_boundary": "absorbing at target",
        "analytic_method": "Appendix 3 finite-sum rearrangement",
        "python_version": platform.python_version(),
        "scipy_version": version("scipy"),
        "matplotlib_version": version("matplotlib")}

    rows = []

    try:
        recurrence = growth_times(
            args.max_target, ratio, args.time_scale)

        for target in range(1, args.max_target + 1):
            analytical = growth_time_analytic(
                target, ratio, args.time_scale)

            difference = abs(recurrence[target] - analytical)

            rows.append({
                "target": target,
                **settings,
                "recurrence_days": recurrence[target],
                "analytic_days": analytical,
                "absolute_difference_days": difference,
                "relative_difference": (
                    difference / analytical
                    if analytical > 0 else None),
                "log10_recurrence_days": (
                    math.log10(recurrence[target])
                    if recurrence[target] > 0 else None),
                "log10_analytic_days": (
                    math.log10(analytical)
                    if analytical > 0 else None)})

    except ValueError as error:
        parser.error(str(error))

    # Export all targets, including the zero-time result at target 1.
    prefix = args.output_prefix

    with open(
        f"{prefix}.csv", "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print("\nFigure A1 comparison:")
    print(
        f"Starting population: 1; "
        f"target abundances: 1 to {args.max_target}")
    print("Lower boundary: m_0=0, equivalent holding at 1")
    print(f"Death-to-birth ratio (q/p): {ratio}")
    print(f"Birth probability: {p}; death probability: {q}")
    print(f"Time constant delta_t: {args.time_scale} days")

    max_difference = max(
        row["absolute_difference_days"] for row in rows)

    print(
        f"Maximum absolute difference: "
        f"{max_difference:.3e} days")
    print(f"Saved {len(rows)} data rows: {prefix}.csv")

    # Target 1 remains in the CSV but has no finite logarithmic value.
    plotted = rows[1:]

    targets = [row["target"] for row in plotted]
    log_analytic = [
        row["log10_analytic_days"] for row in plotted]
    log_recurrence = [
        row["log10_recurrence_days"] for row in plotted]

    fig, ax = plt.subplots(
        figsize=(3.7, 3.5), layout="constrained")

    point_marker = "o" if len(plotted) == 1 else None

    ax.plot(
        targets, log_analytic,
        color="#C43C39",
        linewidth=1.8,
        marker=point_marker,
        markersize=3,
        label="Analytical",
        zorder=3)

    ax.plot(
        targets, log_recurrence,
        color="black",
        linewidth=1.3,
        linestyle=(0, (5, 3)),
        marker=point_marker,
        markersize=2,
        label="Recurrence",
        zorder=4)

    ax.set_xlabel(
        r"$n$ Target population size", labelpad=7)
    ax.set_ylabel(
        r"$\log_{10}(M_{1,n}/\mathrm{day})$", labelpad=8)

    ax.set_xlim(0, args.max_target)

    all_logs = log_analytic + log_recurrence
    padding = max(
        0.08, 0.05 * (max(all_logs) - min(all_logs)))

    ax.set_ylim(
        min(all_logs) - padding, max(all_logs) + padding)

    ax.xaxis.set_major_locator(
        MaxNLocator(nbins=5, integer=True))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.set_facecolor("#FCFCFC")

    for spine in ax.spines.values():
        spine.set_color("#777777")
        spine.set_linewidth(0.7)

    ax.tick_params(
        which="major", direction="out", length=3, width=0.7)
    ax.tick_params(which="minor", length=0)

    ax.set_axisbelow(True)
    ax.grid(
        which="major", color="#AFAFAF", linewidth=0.65)
    ax.grid(
        which="minor", color="#D6D6D6", linewidth=0.45)

    ax.legend(
        loc="upper left", frameon=False, fontsize=8)

    ax.text(
        0.97, 0.07,
        r"$i=1$" + "\n"
        + rf"$\frac{{q}}{{p}}={ratio:.8f}$" + "\n"
        + rf"$\Delta t={args.time_scale:.6g}\ \mathrm{{day}}$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        bbox={
            "facecolor": "#FCFCFC",
            "edgecolor": "none",
            "pad": 3})

    for extension in ["png", "pdf", "svg"]:
        fig.savefig(
            f"{prefix}.{extension}",
            dpi=600,
            facecolor="white")

    print(
        f"Saved plots: {prefix}.png, "
        f"{prefix}.pdf, {prefix}.svg")

    plt.show()    

#Figure A2: comparison of extinction-duration calculations

if args.figure == "A2":
    if args.max_initial < 1:
        parser.error("--max-initial must be at least 1.")

    # Parameters and provenance retained in every CSV row.
    settings = {
        "ratio_q_over_p": ratio,
        "p": p,
        "q": q,
        "delta_t_days": args.time_scale,
        "stopping_abundance": 0,
        "upper_boundary": "unbounded",
        "analytic_method": "Appendix 4 equivalent integral",
        "python_version": platform.python_version(),
        "scipy_version": version("scipy"),
        "matplotlib_version": version("matplotlib")}

    rows = []

    try:
        recurrence = extinction_times(
            args.max_initial, ratio, args.time_scale)

        for starting in range(1, args.max_initial + 1):
            analytical, error = extinction_time_analytic(
                starting, ratio, args.time_scale)

            difference = abs(recurrence[starting] - analytical)

            rows.append({
                "initial": starting,
                **settings,
                "recurrence_days": recurrence[starting],
                "analytic_days": analytical,
                "quadrature_error_days": error,
                "absolute_difference_days": difference,
                "relative_difference": difference / analytical,
                "log10_recurrence_days": math.log10(
                    recurrence[starting]
                ),
                "log10_analytic_days": math.log10(analytical)})

    except ValueError as error:
        parser.error(str(error))

    # Export the same values used in the plot.
    prefix = args.output_prefix

    with open(
        f"{prefix}.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print("\nFigure A2 comparison:")
    print(f"Starting abundances: 1 to {args.max_initial}")
    print(f"Death-to-birth ratio (q/p): {ratio}")
    print(f"Birth probability: {p}; death probability: {q}")
    print(f"Time constant delta_t: {args.time_scale} days")
    print("Stopping abundance: 0; upper boundary: none")

    max_difference = max(
        row["absolute_difference_days"] for row in rows)

    print(
        f"Maximum absolute difference: "
        f"{max_difference:.3e} days")
    print(f"Saved {len(rows)} data rows: {prefix}.csv")

    # Extract the plotting coordinates.
    abundances = [row["initial"] for row in rows]
    log_analytic = [row["log10_analytic_days"] for row in rows]
    log_recurrence = [
        row["log10_recurrence_days"] for row in rows]

    fig, ax = plt.subplots(
        figsize=(3.7, 3.5),
        layout="constrained")

    # A marker also makes a one-abundance plot visible.
    point_marker = "o" if args.max_initial == 1 else None

    ax.plot(
        abundances, log_analytic,
        color="#C43C39",
        linewidth=1.8,
        marker=point_marker,
        markersize=3,
        label="Analytical",
        zorder=3)

    ax.plot(
        abundances, log_recurrence,
        color="black",
        linewidth=1.3,
        linestyle=(0, (5, 3)),
        marker=point_marker,
        markersize=2,
        label="Recurrence",
        zorder=4)

    ax.set_xlabel(
        r"$n$ Initial population size",
        labelpad=7)
    ax.set_ylabel(
        r"$\log_{10}(T_n/\mathrm{day})$",
        labelpad=8)

    ax.set_xlim(0, args.max_initial)

    # Adapt the vertical limits to the selected parameters.
    all_logs = log_analytic + log_recurrence
    padding = max(
        0.08,
        0.05 * (max(all_logs) - min(all_logs)))

    ax.set_ylim(
        min(all_logs) - padding,
        max(all_logs) + padding)

    ax.xaxis.set_major_locator(
        MaxNLocator(nbins=5, integer=True))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.set_facecolor("#FCFCFC")

    for spine in ax.spines.values():
        spine.set_color("#777777")
        spine.set_linewidth(0.7)

    ax.tick_params(
        which="major",
        direction="out",
        length=3,
        width=0.7)
    ax.tick_params(which="minor", length=0)

    ax.set_axisbelow(True)
    ax.grid(
        which="major",
        color="#AFAFAF",
        linewidth=0.65)
    ax.grid(
        which="minor",
        color="#D6D6D6",
        linewidth=0.45)

    ax.legend(
        loc="upper left",
        frameon=False,
        fontsize=8)

    ax.text(
        0.97, 0.07,
        rf"$\frac{{q}}{{p}}={ratio:.8f}$" + "\n"
        + rf"$\Delta t={args.time_scale:.6g}\ \mathrm{{day}}$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        bbox={
            "facecolor": "#FCFCFC",
            "edgecolor": "none",
            "pad": 3})

    for extension in ["png", "pdf", "svg"]:
        fig.savefig(
            f"{prefix}.{extension}",
            dpi=600,
            facecolor="white")

    print(
        f"Saved plots: {prefix}.png, "
        f"{prefix}.pdf, {prefix}.svg")

    plt.show()

# Figure A3: caption-consistent parameter noise, fixed per curve.

if args.figure == "A3":
    if args.max_initial < 1 or args.draws < 2 or args.show_draws < 1:
        parser.error(
            "Require max-initial >= 1, draws >= 2 and show-draws >= 1.")

    if args.seed < 0:
        parser.error("The seed must be non-negative.")

    if not math.isfinite(args.noise_halfwidth) or args.noise_halfwidth < 0:
        parser.error(
            "Noise half-width must be finite and non-negative.")

    base_q = q
    width = args.noise_halfwidth

    if not 0.5 < base_q < 1:
        parser.error(
            "A3 requires a baseline death probability between 0.5 and 1.")

    if base_q - width < 0.5 or base_q + width >= 1:
        parser.error(
            "This finite-mean mode requires the whole noise interval "
            "within [0.5, 1).")

    # Draw once per curve, then hold that probability fixed.
    rng = np.random.default_rng(args.seed)

    sampled_q = base_q + rng.uniform(
        -width, width, size=args.draws)

    if np.any(sampled_q <= 0.5):
        parser.error(
            "A sampled q <= 0.5 has no finite unbounded mean. "
            "No values were redrawn.")

    # Derive ratios consistently for the baseline and parameter draws.
    q_values = [base_q] + [float(value) for value in sampled_q]
    curves = []

    try:
        for death in q_values:
            current_ratio = death / (1 - death)

            curves.append(extinction_times(
                args.max_initial, current_ratio, args.time_scale))

    except ValueError as error:
        parser.error(str(error))

    # Exclude T_0 from the arrays used for plotting and summaries.
    baseline = np.asarray(curves[0][1:])
    noisy = np.asarray([curve[1:] for curve in curves[1:]])

    # Calculate statistics in days, before taking logarithms.
    mean_days = noisy.mean(axis=0)
    sd_days = noisy.std(axis=0, ddof=1)

    low_days, median_days, high_days = np.quantile(
        noisy, [0.025, 0.5, 0.975],
        axis=0, method="linear")

    shown = min(args.show_draws, args.draws)
    prefix = args.output_prefix

    settings = {
        "baseline_q": base_q,
        "requested_baseline_ratio": ratio,
        "delta_t_days": args.time_scale,
        "noise_halfwidth_q": width,
        "noise_draws": args.draws,
        "shown_draws": shown,
        "seed": args.seed,
        "rng_algorithm": type(rng.bit_generator).__name__,
        "noise_rule": "one fixed q draw per curve; p=1-q",
        "initialisation": "perturbed p used consistently",
        "boundary": "0 absorbing; no upper limit",
        "python_version": platform.python_version(),
        "numpy_version": version("numpy"),
        "scipy_version": version("scipy"),
        "matplotlib_version": version("matplotlib")}

    # Export every curve, including the baseline as draw_id 0.
    if not args.summary_only:
        fields = [
            "draw_id", "kind", "initial", "p", "q", "ratio_q_over_p",
            "delta_q", "mean_extinction_days",
            "log10_mean_extinction_days"] + list(settings)

        with open(
            f"{prefix}.csv", "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()

            for draw_id, (death, curve) in enumerate(zip(q_values, curves)):
                for initial in range(1, args.max_initial + 1):
                    writer.writerow({
                        "draw_id": draw_id,
                        "kind": (
                            "baseline" if draw_id == 0 else "parameter_draw"),
                        "initial": initial,
                        "p": 1 - death,
                        "q": death,
                        "ratio_q_over_p": death / (1 - death),
                        "delta_q": death - base_q,
                        "mean_extinction_days": curve[initial],
                        "log10_mean_extinction_days": math.log10(
                            curve[initial]),
                        **settings})

    # Summarise variation across parameter draws, excluding the baseline.
    fields = [
        "initial", "baseline_days", "noise_mean_days",
        "noise_median_days", "noise_sample_sd_days",
        "p025_days", "p975_days",
        "mean_minus_baseline_days"] + list(settings)

    with open(
        f"{prefix}_summary.csv", "w",
        newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        for index in range(args.max_initial):
            writer.writerow({
                "initial": index + 1,
                "baseline_days": baseline[index],
                "noise_mean_days": mean_days[index],
                "noise_median_days": median_days[index],
                "noise_sample_sd_days": sd_days[index],
                "p025_days": low_days[index],
                "p975_days": high_days[index],
                "mean_minus_baseline_days": (
                    mean_days[index] - baseline[index]),
                **settings})

    print("\nFigure A3: fixed parameter draws")
    print(f"Baseline q: {base_q}; baseline p: {1 - base_q}")
    print(f"Uniform noise half-width in q: {width}")
    print(f"Draws: {args.draws}; shown: {shown}; seed: {args.seed}")
    print(
        f"Time constant: {args.time_scale} days; "
        f"n = 1 to {args.max_initial}")
    print(
        f"At n={args.max_initial}: "
        f"baseline mean = {baseline[-1]:.6f} days")
    print(f"Across parameter draws: mean = {mean_days[-1]:.6f} days")
    print(f"Sample SD = {sd_days[-1]:.6f} days")
    print(
        "Empirical 2.5–97.5 percentile range = "
        f"{low_days[-1]:.6f} to {high_days[-1]:.6f} days")

    abundances = np.arange(1, args.max_initial + 1)
    log_noisy = np.log10(noisy)
    log_baseline = np.log10(baseline)

    fig, ax = plt.subplots(
        figsize=(3.7, 3.5), layout="constrained")

    marker = "o" if args.max_initial == 1 else None

    for index in range(shown):
        ax.plot(
            abundances, log_noisy[index],
            color="#D1495B",
            linewidth=1.0,
            alpha=0.25,
            marker=marker,
            markersize=2,
            zorder=2)

    ax.plot(
        abundances, log_baseline,
        color="black",
        linewidth=1.8,
        marker=marker,
        markersize=3,
        label="Baseline",
        zorder=4)

    # Legend key only: empty arrays add no plotted data.
    ax.plot(
        [], [], color="#D1495B",
        linewidth=1.0, label="Parameter draws")

    ax.set_xlabel(
        r"$n$ Initial population size", labelpad=7)
    ax.set_ylabel(
        r"$\log_{10}(T_n/\mathrm{day})$", labelpad=8)

    ax.set_xlim(0, args.max_initial)

    bottom = min(
        float(log_noisy[:shown].min()),
        float(log_baseline.min()))
    top = max(
        float(log_noisy[:shown].max()),
        float(log_baseline.max()))

    padding = max(0.08, 0.05 * (top - bottom))

    ax.set_ylim(
        0.5 * math.floor((bottom - padding) / 0.5),
        0.5 * math.ceil((top + padding) / 0.5))

    ax.xaxis.set_major_locator(
        MaxNLocator(nbins=5, integer=True))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(
        MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.set_facecolor("#FCFCFC")

    for spine in ax.spines.values():
        spine.set_color("#777777")
        spine.set_linewidth(0.7)

    ax.tick_params(
        which="major", direction="out", length=3, width=0.7)
    ax.tick_params(which="minor", length=0)
    ax.set_axisbelow(True)

    ax.grid(
        which="major", color="#AFAFAF", linewidth=0.65)
    ax.grid(
        which="minor", color="#D6D6D6", linewidth=0.45)

    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0, 1.02),
        borderaxespad=0,
        ncol=2,
        frameon=False,
        fontsize=7.5,
        columnspacing=1.0,
        handlelength=1.8)

    ax.text(
        0.97, 0.07,
        rf"$q_0={base_q:.6f}$" + "\n"
        + rf"$\epsilon\sim U(-{width:.4g}, {width:.4g})$" + "\n"
        + rf"$\Delta t={args.time_scale:.6g}\ \mathrm{{day}}$" + "\n"
        + f"{shown}/{args.draws} draws shown; seed={args.seed}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        bbox={
            "facecolor": "#FCFCFC",
            "edgecolor": "none",
            "pad": 3})

    for extension in ["png", "pdf", "svg"]:
        fig.savefig(
            f"{prefix}.{extension}",
            dpi=600,
            facecolor="white")

    if not args.summary_only:
        print(f"Saved full curve data: {prefix}.csv")

    print(
        f"Saved {prefix}_summary.csv "
        "and PNG/PDF/SVG plots.")

    plt.show()