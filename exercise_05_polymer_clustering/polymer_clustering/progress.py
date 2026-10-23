# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 21 October 2026

"""Observe execution without adding model events, observations, or random draws.

Terminal rendering and the small worker status file are deliberately outside the
scientific settings. A progress-file or terminal failure must not lose a run.
"""

import json
from pathlib import Path
import shutil
import sys
import time


def add_progress_argument(parser) -> None:
    """Share one display option between single-run and batch entry points."""
    parser.add_argument("--progress", choices=("auto", "on", "off"), default="auto",
                        help="Progress: auto in a terminal, on also in logs, or off. "
                             "Display only; does not change saved model settings.")


def bar(fraction: float) -> str:
    """Use ASCII characters so bars work on ordinary Linux terminals."""
    fraction = min(1.0, max(0.0, fraction))
    filled = int(12 * fraction)
    return "[" + "#" * filled + "-" * (12 - filled) + "]"


class ConsoleProgress:
    """Redraw one TTY line; explicit non-TTY output uses occasional plain lines."""

    def __init__(self, mode="auto", stream=None, clock=None):
        self.stream = sys.stderr if stream is None else stream
        self.clock = time.monotonic if clock is None else clock
        self.tty = bool(getattr(self.stream, "isatty", lambda: False)())
        self.enabled = mode == "on" or (mode == "auto" and self.tty)
        self.interval = 1.0 if self.tty else 10.0
        self.last_draw = float("-inf")
        self.line_width = 0

    def _write(self, value: str) -> None:
        try:
            self.stream.write(value)
            self.stream.flush()
        except (OSError, ValueError):
            # For example, a closed pipe should disable progress, not the SSA.
            self.enabled = False
            self.line_width = 0

    def update(self, value: str, force=False) -> None:
        if not self.enabled:
            return
        now = self.clock()
        if not force and now - self.last_draw < self.interval:
            return
        self.last_draw = now
        if self.tty:
            # Avoid wrapping: overwriting a wrapped line would leave stale text.
            width = max(1, shutil.get_terminal_size((100, 24)).columns - 1)
            value = value[:width]
            padding = " " * max(0, min(width, self.line_width) - len(value))
            self._write("\r" + value + padding)
            self.line_width = len(value) if self.enabled else 0
        else:
            self._write(value + "\n")

    def close(self) -> None:
        if self.line_width:
            self._write("\n")
            self.line_width = 0

    def message(self, value: str) -> None:
        """End the live line before printing a normal lifecycle message."""
        self.close()
        print(value, file=self.stream, flush=True)


class RunProgress:
    """Throttle read-only model-time callbacks and publish atomic worker status."""

    def __init__(self, duration: float, mode="auto", status_file: Path | None = None,
                 stream=None, clock=None):
        self.console = ConsoleProgress(mode, stream=stream, clock=clock)
        self.clock = self.console.clock
        self.started = self.clock()
        self.last_publish = float("-inf")
        self.status_file = status_file
        self.duration = duration
        self.model_time = 0.0
        self.events = 0
        self.phase_name = "initialising"

    @property
    def callback(self):
        # With display off and no parent status file, the solver has no callback.
        return self.observe if self.console.enabled or self.status_file is not None else None

    def phase(self, name: str) -> None:
        self.phase_name = name
        self.publish(force=True)

    def observe(self, model_time: float, events: int) -> None:
        self.model_time, self.events = model_time, events
        self.publish()

    def publish(self, force=False) -> None:
        if not self.console.enabled and self.status_file is None:
            return
        now = self.clock()
        # Workers update at most once per real second, even with many SSA events.
        if not force and now - self.last_publish < 1.0:
            return
        self.last_publish = now
        elapsed = now - self.started
        fraction = min(1.0, max(0.0, self.model_time / self.duration)) if self.duration else 0.0
        record = dict(schema_version=1, phase=self.phase_name, model_time=self.model_time,
                      duration=self.duration, events=self.events, wall_seconds=elapsed)
        if self.status_file is not None:
            try:
                temporary = self.status_file.with_suffix(self.status_file.suffix + ".tmp")
                temporary.write_text(json.dumps(record, allow_nan=False) + "\n", encoding="utf-8")
                temporary.replace(self.status_file)
            except OSError as error:
                self.status_file = None
                self.console.close()
                # Reporting remains best-effort even if stderr is unavailable.
                self.console._write(f"Progress file unavailable: {error}\n")
        text = (f"{self.phase_name}: {bar(fraction)} simulation {fraction:.1%}; "
                f"t={self.model_time:.6g}/{self.duration:g} a.u.; "
                f"events={self.events:,}; elapsed={elapsed:.0f}s") if self.duration else (
                    f"{self.phase_name}; elapsed={elapsed:.0f}s")
        self.console.update(text, force=force)

    def finish(self, status: str) -> None:
        self.phase(status)
        self.console.close()


class BatchProgress:
    """Count verified complete jobs and show each running child's latest phase.

The bar measures jobs, not a runtime forecast. Simulated time is not proportional
to CPU work, and a worker may still be writing GIFs after its simulation ends.
"""

    def __init__(self, directory: Path, jobs: list[dict], mode="auto"):
        self.directory, self.jobs = directory, jobs
        self.console = ConsoleProgress(mode, stream=sys.stdout)
        self.started = time.monotonic()

    def update(self, force=False) -> None:
        if not self.console.enabled:
            return
        # Do not read worker files on every 0.2-second process poll.
        now = time.monotonic()
        if not force and now - self.console.last_draw < self.console.interval:
            return
        complete = sum(j["status"] == "complete" for j in self.jobs)
        total = len(self.jobs)
        parts = [f"Batch {bar(complete / total)} {complete}/{total} complete"]
        for job in self.jobs:
            if job["status"] != "running":
                continue
            phase = "starting"
            try:
                attempt = self.directory / job["attempts"][-1]["directory"]
                record = json.loads((attempt / "progress.json").read_text(encoding="utf-8"))
                phase = str(record["phase"])
                if phase == "simulating":
                    fraction = min(1.0, max(0.0, record["model_time"] / record["duration"]))
                    phase = f"simulation {fraction:.1%}"
            except (OSError, ValueError, KeyError, TypeError, ZeroDivisionError):
                # A child can be importing, or its optional status can be absent.
                pass
            parts.append(f"{job['id']}: {phase}")
        for status in ("incomplete", "failed", "interrupted"):
            count = sum(j["status"] == status for j in self.jobs)
            if count:
                parts.append(f"{count} {status}")
        # Put elapsed time last so a normal 80-column terminal shows at least
        # the first active job's percentage before truncating secondary details.
        parts.append(f"elapsed={now - self.started:.0f}s")
        self.console.update(" | ".join(parts), force=force)

    def message(self, value: str) -> None:
        self.console.message(value)

    def close(self) -> None:
        self.update(force=True)
        self.console.close()
