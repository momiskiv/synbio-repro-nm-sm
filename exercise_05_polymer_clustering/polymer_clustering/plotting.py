# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""One scoped typography policy for every static plot and animation frame."""

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from functools import wraps
import hashlib
from pathlib import Path
from typing import Callable, Iterator, ParamSpec, TypeVar
import warnings

import matplotlib
from matplotlib import font_manager


@dataclass(frozen=True)
class PlotConfig:
    """Presentation settings are saved for replay, separately from the biology."""

    font_family: str = "Arial"
    font_path: str | None = None
    strict_font: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.font_family, str) or not self.font_family.strip():
            raise ValueError("font_family must be a non-empty font family name.")
        if self.font_path is not None and (
            not isinstance(self.font_path, str) or not self.font_path.strip()
        ):
            raise ValueError("font_path must be a non-empty file path or null.")
        if type(self.strict_font) is not bool:
            raise ValueError("strict_font must be true or false.")


@dataclass(frozen=True)
class ResolvedFont:
    """Record the font actually selected, including any clearly identified fallback."""

    requested: str
    family: str
    path: Path
    fallback_used: bool
    registered_path: Path | None = None

    def metadata(self) -> dict[str, object]:
        return {
            "requested_font_family": self.requested,
            "resolved_font_family": self.family,
            "font_file": str(self.path),
            "font_file_sha256": hashlib.sha256(self.path.read_bytes()).hexdigest(),
            "fallback_used": self.fallback_used,
            "registered_font_file": str(self.registered_path) if self.registered_path else None,
            "svg_font_type": "path",
        }

    def notice(self) -> str:
        return (f"Font '{self.requested}' was not found; using '{self.family}'. "
                "Install the requested font or supply --font-path; "
                "--strict-font requires the requested font to be available.")

    def rc_params(self) -> dict[str, object]:
        # A single resolved family avoids repeated missing-font log messages.
        # SVG glyph paths retain their appearance when opened on another machine.
        return {
            "font.family": [self.family], "font.sans-serif": [self.family],
            "mathtext.fontset": "custom", "mathtext.rm": self.family,
            "mathtext.sf": self.family, "mathtext.tt": self.family,
            "mathtext.cal": self.family, "mathtext.it": f"{self.family}:italic",
            "mathtext.bf": f"{self.family}:bold", "svg.fonttype": "path",
        }


def resolve_font(config: PlotConfig) -> ResolvedFont:
    """Use installed Arial by default, or register a user-supplied font file.

    A supplied file must match the requested family. Arial is not bundled or
    downloaded by the project. If unavailable, an announced Liberation Sans or
    bundled DejaVu Sans fallback keeps simulations usable unless strict mode is on.
    """
    registered = None
    if config.font_path is not None:
        registered = Path(config.font_path).expanduser().resolve()
        if not registered.is_file():
            raise ValueError(f"Font file does not exist: {registered}")
        try:
            family = font_manager.FontProperties(fname=registered).get_name()
            if family.casefold() != config.font_family.strip().casefold():
                raise ValueError(
                    f"Font file family is '{family}', but font_family is '{config.font_family}'."
                )
            font_manager.fontManager.addfont(registered)
        except (RuntimeError, OSError) as error:
            raise ValueError(f"Could not load font file: {registered}: {error}") from error
    requested = config.font_family.strip()
    try:
        path = Path(font_manager.findfont(font_manager.FontProperties(family=[requested]),
                                          fallback_to_default=False))
        fallback = False
    except ValueError as error:
        if config.strict_font:
            raise ValueError(
                f"Requested font '{requested}' is unavailable; install it or use --font-path."
            ) from error
        for family in ("Liberation Sans", "DejaVu Sans"):
            try:
                path = Path(font_manager.findfont(font_manager.FontProperties(family=[family]),
                                                  fallback_to_default=False))
                break
            except ValueError:
                continue
        else:
            raise ValueError("No supported fallback font was found.") from error
        fallback = True
    family = font_manager.FontProperties(fname=path).get_name()
    return ResolvedFont(requested, family, path, fallback, registered)


_ACTIVE_FONT: ContextVar[ResolvedFont | None] = ContextVar("polymer_plot_font", default=None)


@contextmanager
def plot_style(font: ResolvedFont | None = None) -> Iterator[ResolvedFont]:
    """Apply typography locally, restoring callers' Matplotlib settings on exit."""
    active = _ACTIVE_FONT.get()
    if font is None and active is not None:
        yield active
        return
    if font is None:
        font = resolve_font(PlotConfig())
        if font.fallback_used:
            warnings.warn(font.notice(), RuntimeWarning)
    token = _ACTIVE_FONT.set(font)
    try:
        with matplotlib.rc_context(font.rc_params()):
            yield font
    finally:
        _ACTIVE_FONT.reset(token)


P = ParamSpec("P")
R = TypeVar("R")


def styled_plot(function: Callable[P, R]) -> Callable[P, R]:
    """Give direct library plot calls the same default policy as CLI-generated plots."""
    @wraps(function)
    def wrapped(*args: P.args, **kwargs: P.kwargs) -> R:
        with plot_style():
            return function(*args, **kwargs)
    return wrapped
