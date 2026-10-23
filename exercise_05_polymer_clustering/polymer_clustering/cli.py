# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 22 October 2026

"""Command-line controls for staged diagnostics and the coupled population model."""

import argparse
import json
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .binding import BindingConfig, CROWDING_THRESHOLD, simulate_binding
from .binding_output import save_binding
from .model import InitialConfig, initialise
from .movement import MovementConfig, simulate_movement
from .movement_output import save_movement
from .output import create_run_directory, save_plot, save_state, software_metadata, write_json
from .plotting import PlotConfig, plot_style, resolve_font
from .progress import RunProgress, add_progress_argument
from .timing import RunTimer, timing_summary
from .population import PopulationConfig
from .population_output import save_population
from .response import ResponseConfig
from .response_output import save_response


def parse_ps_threshold(value: str) -> int | None:
    """Distinguish a positive threshold from an explicitly disabled PS cutoff."""
    if value.casefold() in ("off", "none"):
        return None
    try:
        threshold = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("Use a positive integer or 'off'.") from error
    if threshold < 1:
        raise argparse.ArgumentTypeError("Use a positive integer or 'off'.")
    return threshold


# User controls: argparse handles --help and converts option values to types.
def build_parser() -> argparse.ArgumentParser:
    """Describe controls without requiring users to edit Python source."""
    parser = argparse.ArgumentParser(
        description=(
            "Create an initial state, run fixed-position binding, or couple binding "
            "and movement, or add a downstream expression/light response. "
            "Full mode adds age-dependent division and strain-specific signal production."
        ),
        allow_abbrev=False,
    )
    # Model options default to None here so an omitted command-line option
    # cannot accidentally overwrite a setting supplied in the JSON file.
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--config", type=Path, help="Read settings from a JSON file.")
    parser.add_argument("--mode", choices=("initial", "binding", "movement", "response", "full"),
                        help="Run a diagnostic stage, or full mode with response, growth and production.")
    parser.add_argument("--grid-size", type=int, help="Square grid side length (default: 50).")
    parser.add_argument("--bacteria", type=int, help="Initial bacteria, 0 to 80 (default: 25).")
    parser.add_argument("--polymers", type=int, help="Initial free polymer count (default: 0).")
    parser.add_argument("--strain", choices=("MM32", "BB170"), help="Strain (default: MM32).")
    parser.add_argument("--seed", type=int, help="Random seed; default: 2026 (our choice).")
    parser.add_argument("--rng-scheme", choices=("independent_v1", "legacy_v02"),
                        help="Placement RNG streams (default: independent_v1).")
    parser.add_argument("--duration", type=float, help="Binding duration in a.u. (default: 100).")
    parser.add_argument("--sample-interval", type=float,
                        help="Binding observation interval in a.u. (default: 5).")
    parser.add_argument("--max-events", type=int,
                        help="Stop and mark incomplete at this event cap (default: 500000).")
    for reaction in ("bs", "pb", "ps"):
        for direction in ("on", "off"):
            parser.add_argument(f"--k-{reaction}-{direction}", type=float,
                                help=f"{reaction.upper()} {direction} rate; see model notes.")
    parser.add_argument("--ps-crowding-threshold", type=parse_ps_threshold,
                        default=argparse.SUPPRESS,
                        help="Block new PS binding at this local polymer count; default 10. "
                             "Use 'off' for a labelled sensitivity intervention.")
    for flag in ("d-s", "d-p", "d-b", "obstruction-factor", "crowding-lambda"):
        parser.add_argument(f"--{flag}", type=float, help="Movement control; see model notes.")
    parser.add_argument("--advection", action=argparse.BooleanOptionalAction, default=None,
                        help="Enable/disable cluster displacement of free chemicals (default: on).")
    parser.add_argument("--validate-every", type=int,
                        help="Full state check after this many events (default: 1000).")
    parser.add_argument("--spatial-frames", type=int,
                        help="Number of spatial animation frames, 2-201 (default: 21).")
    parser.add_argument("--event-log", choices=("compact", "all"),
                        help="Movement log: compact omits signal hops; all logs every event.")
    for flag in ("k-syn", "k-deg", "response-initial", "response-midpoint", "light-scale",
                 "light-map-max"):
        parser.add_argument(f"--{flag}", type=float, help="Response control; see model notes.")
    parser.add_argument("--growth", action=argparse.BooleanOptionalAction, default=None,
                        help="Full mode: enable/disable division (default: enabled).")
    for flag in ("doubling-time", "division-width", "division-retry", "production-rate",
                 "daughter-expression"):
        parser.add_argument(f"--{flag}", type=float, help="Population control; see model notes.")
    parser.add_argument("--population-cap", type=int, help="Maximum bacteria, 0-80 (default: 80).")
    parser.add_argument("--font-family", help="Plot font family (default: Arial).")
    parser.add_argument("--font-path", help="Optional local font file matching --font-family.")
    parser.add_argument("--strict-font", action=argparse.BooleanOptionalAction, default=None,
                        help="Require the requested font instead of an announced fallback.")
    parser.add_argument(
        "--output", type=Path, default=Path("outputs"),
        help="Parent output folder, relative to your working directory (default: outputs).",
    )
    add_progress_argument(parser)
    # Private worker-to-controller status, kept outside scientific parameters.
    parser.add_argument("--progress-file", type=Path, help=argparse.SUPPRESS)
    return parser


# Run lifecycle: resolve settings, create outputs, and record the outcome.
def main(argv: list[str] | None = None) -> int:
    """Validate settings, generate a state, and preserve success or failure metadata."""
    parser = build_parser()
    args = parser.parse_args(argv)
    # Elapsed real time is independent of model time and UTC clock adjustments.
    timer = RunTimer()
    started_utc = datetime.now(timezone.utc).isoformat()
    # Precedence is defaults < JSON configuration < explicit command-line flags.
    # asdict converts the configuration dataclass to a plain dictionary.
    initial_keys = set(asdict(InitialConfig()))
    binding_keys = set(asdict(BindingConfig()))
    movement_keys = set(asdict(MovementConfig()))
    response_keys = set(asdict(ResponseConfig()))
    population_keys = set(asdict(PopulationConfig()))
    plotting_keys = set(asdict(PlotConfig()))
    settings = {**asdict(InitialConfig()), **asdict(BindingConfig()),
                **asdict(MovementConfig()), **asdict(ResponseConfig()),
                **asdict(PopulationConfig()), **asdict(PlotConfig()), "mode": "initial"}
    try:
        if args.config is not None:
            supplied = json.loads(args.config.read_text(encoding="utf-8"))
            if not isinstance(supplied, dict):
                raise ValueError("The configuration must be a JSON object.")
            # Reject unexpected keys so a mistyped setting cannot be ignored.
            unknown = set(supplied) - set(settings)
            if unknown:
                raise ValueError(f"Unknown configuration keys: {', '.join(sorted(unknown))}")
            settings.update(supplied)
            if "rng_scheme" not in supplied:
                # Old settings files must retain the old placement draw order.
                settings["rng_scheme"] = "legacy_v02"
        for key in settings:
            override = getattr(args, key, None)
            if override is not None:
                settings[key] = override
        # An explicit '--ps-crowding-threshold off' must override JSON even
        # though its value is None. Omission leaves defaults/JSON untouched.
        if hasattr(args, "ps_crowding_threshold"):
            settings["ps_crowding_threshold"] = args.ps_crowding_threshold
        # Unpack the merged dictionary as named arguments, then validate it.
        config = InitialConfig(**{key: settings[key] for key in initial_keys})
        binding = BindingConfig(**{key: settings[key] for key in binding_keys})
        movement = MovementConfig(**{key: settings[key] for key in movement_keys})
        response = ResponseConfig(**{key: settings[key] for key in response_keys})
        population = PopulationConfig(**{key: settings[key] for key in population_keys})
        plotting = PlotConfig(**{key: settings[key] for key in plotting_keys})
        # Invalid font files/strict requests fail before a run folder is created.
        font = resolve_font(plotting)
        if settings["mode"] not in ("initial", "binding", "movement", "response", "full"):
            raise ValueError("mode must be initial, binding, movement, response, or full.")
        full = settings["mode"] == "full"
        responding = settings["mode"] in ("response", "full")
        moving = settings["mode"] in ("movement", "response", "full")
        binding_mode = settings["mode"] != "initial"
        if binding_mode:
            binding.validate_polymer_rates(config.polymers)
        elif any(settings[key] != asdict(BindingConfig())[key] for key in binding_keys):
            raise ValueError("Binding controls require --mode binding.")
        if not moving and any(settings[key] != asdict(MovementConfig())[key]
                              for key in movement_keys):
            raise ValueError("Movement controls require --mode movement.")
        if not responding and any(settings[key] != asdict(ResponseConfig())[key]
                                  for key in response_keys):
            raise ValueError("Response controls require --mode response.")
        if not full and any(settings[key] != asdict(PopulationConfig())[key]
                            for key in population_keys):
            raise ValueError("Population controls require --mode full.")
        if full and config.bacteria > population.population_cap:
            raise ValueError("Initial bacteria exceed population_cap.")
    except (OSError, ValueError, TypeError) as error:
        parser.error(str(error))

    # Invalid settings exit above, before an output directory is created.
    try:
        directory = create_run_directory(args.output, settings["mode"])
    except OSError as error:
        print(f"Could not create the output directory: {error}\n"
              f"{timing_summary(timer.snapshot())}", file=sys.stderr)
        return 1

    if font.fallback_used:
        print(font.notice(), file=sys.stderr)

    # These records explain what was requested and which code performed it.
    # The INIT identifiers refer to the choices in docs/model-notes.md.
    metadata: dict[str, object] = {
        "status": "started",
        "milestone": ("population_coupled_model" if full else
                      "binding_movement_and_response" if responding else
                      "binding_and_movement" if moving else
                      "binding_only_fixed_positions" if binding_mode else "initialisation_only"),
        "started_utc": started_utc,
        "arguments": sys.argv[1:] if argv is None else argv,
        "working_directory": str(Path.cwd()),
        "parameters": ({key: value for key, value in settings.items()
                        if (moving or key not in movement_keys)
                        and (responding or key not in response_keys)
                        and (full or key not in population_keys)}
                       if binding_mode else {**asdict(config), **asdict(plotting)}),
        "initial_total_signals": config.signal_count,
        "source": {"doi": "10.1038/nchem.1793", "pages": ["S20", "S21", "S23"]},
        "interpretations": ["INIT-01", "INIT-02", "INIT-03", "INIT-04"],
        "software": software_metadata(),
        "plotting": font.metadata(),
    }
    metadata["initial_random_streams"] = config.rng_scheme
    if config.rng_scheme == "independent_v1":
        metadata["interpretations"].append("INIT-05")
    if binding_mode:
        metadata["interpretations"] += [f"BIND-0{i}" for i in range(1, 7)]
        metadata["binding_random_stream"] = "PCG64(SeedSequence([seed, 1]))"
        metadata["binding_kinetics"] = binding.rate_metadata()
        metadata["ps_crowding_rule"] = {
            "threshold": binding.ps_crowding_threshold,
            "source_default": CROWDING_THRESHOLD,
            "matches_source_default": binding.ps_crowding_threshold == CROWDING_THRESHOLD,
            "scope": "PS accessibility only; PB binding, transport and ordinary PS release unchanged.",
        }
        if binding.ps_crowding_threshold != CROWDING_THRESHOLD:
            metadata["interpretations"].append("CROWD-01")
        metadata["excluded_processes"] = ["movement", "growth", "signal_production",
                                          "expression", "luminescence"]
    if moving:
        metadata["interpretations"] += [f"MOVE-0{i}" for i in range(1, 8)]
        metadata["source"]["pages"] = ["S19", "S20", "S21", "S22", "S23"]
        metadata["binding_random_stream"] = "PCG64(SeedSequence([seed, 2]))"
        metadata["excluded_processes"] = ["growth", "signal_production",
                                          "expression", "luminescence"]
        metadata["event_log"] = movement.event_log
    if responding:
        metadata["interpretations"] += [f"RESP-0{i}" for i in range(1, 6)]
        metadata["excluded_processes"] = ["growth", "signal_production"]
        metadata["response_model"] = {
            "midpoint": response.midpoint(config.strain), "k_syn": response.k_syn,
            "k_deg": response.k_deg, "initial_expression": response.response_initial,
            "light_scale": response.light_scale,
            "theoretical_light_upper_limit_per_cell": 2 * response.light_scale,
            "decay_rule": "positive degradation magnitude; stop decay at constant sensing",
            "population_light": "raw sum and mean; no baseline subtraction or normalisation",
            "coupling": "downstream only; no feedback into stochastic event rates",
            "light_map_max": response.light_map_max,
        }
    if full:
        metadata["interpretations"] += [f"POP-0{i}" for i in range(1, 7)]
        if population.doubling_time == population.division_width / 2:
            metadata["interpretations"].append("POP-07")
        metadata["excluded_processes"] = []
        metadata["population_model"] = {
            "growth": population.growth,
            "division_distribution": "symmetric triangular lifetime; width is full support",
            "positive_lifetime_policy": "resample finite-precision zero endpoints",
            "division_timing_matches_source_defaults": (population.doubling_time == 15000.0
                                                       and population.division_width == 20000.0),
            "division_support": [population.doubling_time - population.division_width / 2,
                                 population.doubling_time + population.division_width / 2],
            "initial_cell_age": 0,
            "population_cap": population.population_cap,
            "division_retry": population.division_retry,
            "daughter_expression": population.daughter_expression,
            "production_rate_per_cell": population.signal_rate(config.strain),
            "conservation": "free + bacterium-bound + polymer-bound = initial + produced",
            "lifetime_random_stream": "PCG64(SeedSequence([seed, 3]))",
            "daughter_placement_random_stream": "PCG64(SeedSequence([seed, 4]))",
            "source_agreement": "implementation includes unresolved reconstruction choices; not yet a published-figure reproduction",
        }
    total_events = None
    progress = RunProgress(binding.duration if binding_mode else 0.0,
                           args.progress, args.progress_file)
    try:
        progress.phase("initialising")
        # The style context covers every text artist, colourbar, and GIF frame.
        # It restores any surrounding application's Matplotlib defaults on exit.
        with plot_style(font):
            with timer.phase("output"):
                write_json(directory / "run.json", metadata)
                write_json(directory / "parameters.json", metadata["parameters"])
            with timer.phase("initialisation"):
                state = initialise(config)
            with timer.phase("output"):
                progress.phase("saving initial outputs")
                save_state(directory, state)
                save_plot(directory, state, config)
            if binding_mode:
                progress.phase("simulating")
                with timer.phase("simulation"):
                    result = (simulate_movement(state, binding, movement, config.seed,
                                                response=response if responding else None,
                                                strain=config.strain,
                                                population=population if full else None,
                                                progress=progress.callback) if moving
                              else simulate_binding(state, binding, config.seed,
                                                    progress=progress.callback))
                cutoff = binding.ps_crowding_threshold
                # S42's faster BB170 timing is a collection condition. Replace
                # the legacy generic note in the presentation layer, keeping
                # the solver and all numerical arrays unchanged.
                if (full and population.growth and config.strain == "BB170"
                        and population.doubling_time == 10000.0
                        and population.division_width == 20000.0):
                    result.scenario_note = "BB170 faster growth: T=10000 a.u."
                crowding_note = ("" if cutoff == CROWDING_THRESHOLD else
                                 f"Sensitivity: PS crowding cutoff {'off' if cutoff is None else cutoff}")
                if moving and crowding_note:
                    result.scenario_note = "\n".join(filter(None, (result.scenario_note, crowding_note)))
                total_events = result.total_events if moving else len(result.events)
                progress.phase("saving outputs")
                with timer.phase("output"):
                    save_binding(directory, result, moving, scenario_note=crowding_note)
                    if moving:
                        save_movement(directory, result)
                    if responding:
                        save_response(directory, result, response.light_map_max)
                        metadata["response"] = result.response_trajectory[-1]
                    if full:
                        save_population(directory, result)
                        metadata["population"] = result.population_trajectory[-1]
                metadata["binding"] = {"status": result.status, "stop_time": result.time,
                                       "events": total_events,
                                       "recorded_events": len(result.events),
                                       "final_counts": result.trajectory[-1]}
                if moving:
                    metadata["movement_event_counts"] = result.event_counts
            # All requested outputs are written before the final timing snapshot.
            metadata["status"] = "complete"
            if binding_mode and result.status != "complete":
                metadata["status"] = "incomplete"
            metadata["finished_utc"] = datetime.now(timezone.utc).isoformat()
            metadata["timing"] = timer.snapshot(total_events)
            write_json(directory / "run.json", metadata)
        progress.finish(metadata["status"])
    except KeyboardInterrupt:
        metadata.update(status="interrupted", finished_utc=datetime.now(timezone.utc).isoformat(),
                        timing=timer.snapshot(total_events))
        try:
            write_json(directory / "run.json", metadata)
        except OSError:
            pass
        progress.finish("interrupted")
        print(f"Run interrupted. Replay from the original seed.\n"
              f"{timing_summary(metadata['timing'])}\nRun directory: {directory}", file=sys.stderr)
        return 130
    except Exception as error:
        # Each phase records elapsed work in finally, including a failed solver
        # or renderer. Preserve timing and the original error when possible.
        metadata["status"] = "failed"
        metadata["error"] = f"{type(error).__name__}: {error}"
        metadata["finished_utc"] = datetime.now(timezone.utc).isoformat()
        metadata["timing"] = timer.snapshot(total_events)
        try:
            write_json(directory / "run.json", metadata)
        except OSError:
            pass
        progress.finish("failed")
        print(f"Run failed: {error}\n{timing_summary(metadata['timing'])}\n"
              f"Run directory: {directory}", file=sys.stderr)
        return 1
    print(timing_summary(metadata["timing"]))
    if metadata["status"] == "incomplete":
        print(f"Run INCOMPLETE: event limit reached at t={result.time:.6g}. "
              f"Increase --max-events to continue with a fresh replay.\nOutputs: {directory}")
        return 3
    message = ("Binding, movement, division, signal production, and response completed." if full else
               "Binding, movement, expression, and predicted light completed." if responding else
               "Binding and movement completed (no growth or light calculation)." if moving else
               "Fixed-position binding completed (no movement or light calculation)."
               if binding_mode else "Initial state created (no dynamics simulated).")
    print(f"{message}\nOutputs: {directory}")
    # The shell receives 0 for success; failures above return 1 or argparse's 2.
    return 0

