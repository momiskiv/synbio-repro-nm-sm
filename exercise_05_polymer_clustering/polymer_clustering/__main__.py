# Part of Module 1: Synthetic Biology and Bioinfomatics for EngBio4Env
# This code was co-worked using OpenAI Codex by N.Musto
# Last updated: 19 October 2026

"""Support ``python -m polymer_clustering`` from the project directory."""

from .cli import main

# Delegate to the CLI and return its exit status to the calling shell.
if __name__ == "__main__":
    raise SystemExit(main())

